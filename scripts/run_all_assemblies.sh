#!/usr/bin/env bash
# Full-scale launcher: assemble and QC every genotype with the existing
# Snakemake workflow, N genotypes at a time (default 3 x 40 threads = 120).
#
# Each genotype runs as its own Snakemake process (target `per_sample`,
# `--config only_samples=<LINE>`), so a failing genotype never stops or
# invalidates the others, logs are per genotype, and only N genotypes hold
# large intermediates at any time. GNU parallel starts the next genotype as
# soon as one finishes. Rerunning the launcher resumes: completed stages are
# not repeated (Snakemake), completed genotypes finish immediately.
#
# Run from pipeline/. Review first:
#   bash scripts/run_all_assemblies.sh --dry-run
# Then start (survives logout):
#   nohup bash scripts/run_all_assemblies.sh > ../logs/launcher.$(date +%Y%m%d_%H%M%S).out 2>&1 &
# or inside `screen`. Monitor:  cat ../logs/assembly_<stamp>/status.tsv
#                               tail ../logs/assembly_<stamp>/<LINE>.log
#
# Options:
#   --jobs N            genotypes in parallel (default 3)
#   --threads T         threads per genotype (default 40)
#   --mem-gb G          memory budget per genotype for Snakemake scheduling (default 300)
#   --min-free-gb F     refuse to start with less free disk on results/ (default 800)
#   --samples "A B"     only these genotypes (default: all discovered), in this order
#   --resources FILE    resources file (default config/resources/cbsuxm01_3x40.yaml)
#   --dry-run           preflight + `snakemake -n` per genotype; executes nothing
#   --force-disk        start even if free disk is below --min-free-gb (not advised)
set -euo pipefail

JOBS=3; THREADS=40; MEM_GB=300; MIN_FREE_GB=800; SAMPLES=""; DRY=0; FORCE_DISK=0
RESOURCES=config/resources/cbsuxm01_3x40.yaml
while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs) JOBS=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --mem-gb) MEM_GB=$2; shift 2 ;;
        --min-free-gb) MIN_FREE_GB=$2; shift 2 ;;
        --samples) SAMPLES=$2; shift 2 ;;
        --resources) RESOURCES=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        --force-disk) FORCE_DISK=1; shift ;;
        -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")/.."
PIPE=$(pwd)
if command -v module >/dev/null 2>&1; then module load python/3.12.7 >/dev/null 2>&1 || true; fi
unset OMP_NUM_THREADS   # rules that need it set it themselves
STAMP=$(date +%Y%m%d_%H%M%S)
LOGDIR=$(realpath -m "$PIPE/../logs/assembly_$STAMP")
[[ $DRY -eq 1 ]] && LOGDIR="$LOGDIR.dryrun"
mkdir -p "$LOGDIR"
say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOGDIR/launcher.log"; }
die() { say "ERROR: $*"; exit 1; }

# ------------------------------------------------------------------ preflight
say "pipeline: $PIPE (git $(git rev-parse --short HEAD 2>/dev/null || echo '?')$(git diff --quiet 2>/dev/null || echo ', UNCOMMITTED CHANGES'))"
[[ -f "$RESOURCES" ]] || die "resources file not found: $RESOURCES"
command -v snakemake >/dev/null || die "snakemake not found (module load python/3.12.7)"
SMK_VER=$(snakemake --version); [[ ${SMK_VER%%.*} -ge 8 ]] || die "Snakemake >= 8 required (found $SMK_VER)"
command -v parallel >/dev/null || die "GNU parallel not found"
python3 - "$PIPE" <<'EOF' || die "tool check failed (see above)"
import os, shlex, sys, yaml
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
from g2f.config import load_config
with open(os.path.join(sys.argv[1], "config", "config.yaml")) as fh:
    cfg = load_config(yaml.safe_load(fh), sys.argv[1])
bad = []
for name in ("hifiasm", "seqkit", "samtools", "bam2fastq", "pbindex", "meryl", "merqury", "busco"):
    for tok in shlex.split(cfg["tools"][name]["cmd"]):
        if tok.startswith("/") and "," not in tok and not os.path.exists(tok):
            bad.append(f"{name}: missing {tok}")
lin = cfg["local"].get("busco_lineage_dir")
for l in cfg["science"]["qc"]["busco"]["lineages"]:
    if not lin or not os.path.isfile(os.path.join(lin, l, "dataset.cfg")):
        bad.append(f"BUSCO lineage {l} not found under {lin}")
for b in bad:
    print("  " + b, file=sys.stderr)
sys.exit(1 if bad else 0)
EOF

RESDIR=$(python3 -c "import os,yaml;print(os.path.normpath(os.path.join('$PIPE', yaml.safe_load(open('config/config.yaml'))['results_dir'])))")
mkdir -p "$RESDIR"
FREE_GB=$(( $(df --output=avail -k "$RESDIR" | tail -1) / 1024 / 1024 ))
say "free disk on $(df --output=target "$RESDIR" | tail -1 | tr -d ' '): ${FREE_GB} GB (minimum ${MIN_FREE_GB} GB)"
if (( FREE_GB < MIN_FREE_GB )); then
    [[ $FORCE_DISK -eq 1 || $DRY -eq 1 ]] || die "not enough free disk; see docs/full_scale_run_plan.md (or --force-disk)"
    say "WARNING: free disk below minimum"
fi
MEM_AVAIL_GB=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo)
ARC_GB=$(awk '$1=="size"{print int($3/1e9)}' /proc/spl/kstat/zfs/arcstats 2>/dev/null || echo 0)
say "memory: ${MEM_AVAIL_GB} GB available + ${ARC_GB} GB ZFS cache (reclaimable); budget ${JOBS} x ${MEM_GB} GB"
(( JOBS * MEM_GB <= MEM_AVAIL_GB + ARC_GB )) || say "WARNING: memory budget exceeds available memory"
say "load average: $(cut -d' ' -f1-3 /proc/loadavg); threads requested: ${JOBS} x ${THREADS} = $((JOBS * THREADS))"

# Genotypes, largest HiFi input first so the longest assembly does not run last.
if [[ -z "$SAMPLES" ]]; then
    SAMPLES=$(python3 scripts/g2f_tool.py discover --json | python3 -c '
import json, sys
recs = json.load(sys.stdin)
recs = recs["samples"] if isinstance(recs, dict) else recs
size = {r["sample"]: sum(i["size"] for i in r["inputs"] if i["tech"] == "hifi" and i["status"] == "use")
        for r in recs}
print(" ".join(sorted(size, key=lambda s: -size[s])))')
fi
[[ -n "$SAMPLES" ]] || die "no genotypes discovered"
say "genotypes (in start order): $SAMPLES"

COMMON=(--profile profiles/biohpc --rerun-triggers mtime)
say "whole-project dry run (discovery, configuration, DAG) ..."
snakemake per_sample "${COMMON[@]}" --cores "$THREADS" -n --config "resources=$RESOURCES" \
    > "$LOGDIR/dryrun_all.log" 2>&1 \
    || die "dry run failed; see $LOGDIR/dryrun_all.log"
grep -A40 '^Job stats' "$LOGDIR/dryrun_all.log" | tee -a "$LOGDIR/launcher.log" || true

python3 - "$LOGDIR/launch.json" <<EOF
import json, subprocess, sys
git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True).stdout.strip()
json.dump({"started": "$STAMP", "dry_run": bool($DRY), "git_commit": git("rev-parse", "HEAD"),
           "git_dirty": bool(git("status", "--porcelain")), "jobs": $JOBS, "threads": $THREADS,
           "mem_gb_per_genotype": $MEM_GB, "resources": "$RESOURCES", "samples": "$SAMPLES".split(),
           "snakemake": "$SMK_VER", "free_disk_gb": $FREE_GB, "host": "$(hostname)"},
          open(sys.argv[1], "w"), indent=2)
EOF

PER_GENO=(snakemake per_sample "${COMMON[@]}" --cores "$THREADS" --resources "mem_mb=$((MEM_GB * 1000))"
          --keep-going --rerun-incomplete)
if [[ $DRY -eq 1 ]]; then
    for s in $SAMPLES; do
        "${PER_GENO[@]}" --config "resources=$RESOURCES" "only_samples=$s" -n \
            > "$LOGDIR/$s.dryrun.log" 2>&1 || die "dry run for $s failed: $LOGDIR/$s.dryrun.log"
        say "$s: $(grep -m1 -E '^total' "$LOGDIR/$s.dryrun.log" | awk '{print $2}') jobs would run"
    done
    say "DRY RUN complete; nothing was executed. Command per genotype:"
    say "  ${PER_GENO[*]} --config resources=$RESOURCES only_samples=<LINE>"
    exit 0
fi

# ------------------------------------------------------------------ execution
say "creating the shared validation-identity stamp"
snakemake identity_stamp "${COMMON[@]}" --cores 1 --config "resources=$RESOURCES" \
    >> "$LOGDIR/launcher.log" 2>&1 \
    || die "could not create the identity stamp"

printf 'genotype\texit\tlog\n' > "$LOGDIR/status.tsv"
say "starting: ${JOBS} genotypes at a time, ${THREADS} threads each; logs in $LOGDIR"
# shellcheck disable=SC2016
printf '%s\n' $SAMPLES | parallel -j "$JOBS" --halt never --joblog "$LOGDIR/joblog.tsv" \
    'echo "[$(date "+%F %T")] START {}" >> '"$LOGDIR/launcher.log"';
     '"${PER_GENO[*]}"' --config resources='"$RESOURCES"' only_samples={} > '"$LOGDIR"'/{}.log 2>&1;
     rc=$?; printf "%s\t%s\t%s\n" {} $rc '"$LOGDIR"'/{}.log >> '"$LOGDIR"'/status.tsv;
     echo "[$(date "+%F %T")] END {} exit=$rc" >> '"$LOGDIR/launcher.log"'; exit $rc' || true

FAILED=$(awk -F'\t' 'NR>1 && $2!=0{print $1}' "$LOGDIR/status.tsv" | tr '\n' ' ')
if [[ -n "$FAILED" ]]; then
    say "FAILED genotypes: $FAILED (see their logs; the others are complete). Fix, then rerun this"
    say "launcher: completed work is not repeated. Project summary not built."
    exit 1
fi
say "all genotypes completed; building the project summary table"
snakemake all "${COMMON[@]}" --cores 4 --config "resources=$RESOURCES" >> "$LOGDIR/launcher.log" 2>&1 \
    || die "project summary failed"
say "DONE. Summary: $(ls -t "$RESDIR"/summary/assembly_summary.*.tsv | head -1)"
