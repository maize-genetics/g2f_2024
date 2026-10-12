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
#   --threads T         threads per genotype (default 40). Enforced: each running
#                       genotype is pinned (taskset) to its own block of T CPUs,
#                       so the run never uses more than jobs x threads CPUs
#   --mem-gb G          memory budget per genotype for Snakemake scheduling (default 300)
#   --reserve-gb R      free disk always kept in reserve by the disk gate (default 150)
#   --retain-hifi M     prepared HiFi FASTQ retention: keep (default; nothing is
#                       deleted), auto (the disk gate may make it temporary per
#                       genotype), or temp (always delete after use). With auto or
#                       temp, keep the raw BAMs: they become the only local copy.
#   --max-wait-h H      a genotype waiting this long for disk space fails (default 48)
#   --samples "A B"     only these genotypes (default: all discovered), in this order
#   --resources FILE    resources file (default config/resources/cbsuxm01_3x40.yaml)
#   --dry-run           preflight + `snakemake -n` per genotype; executes nothing
set -euo pipefail

JOBS=3; THREADS=40; MEM_GB=300; RESERVE_GB=150; RETAIN_MODE=keep; MAX_WAIT_H=48
SAMPLES=""; DRY=0
RESOURCES=config/resources/cbsuxm01_3x40.yaml
while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs) JOBS=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --mem-gb) MEM_GB=$2; shift 2 ;;
        --reserve-gb) RESERVE_GB=$2; shift 2 ;;
        --retain-hifi) RETAIN_MODE=$2; shift 2 ;;
        --max-wait-h) MAX_WAIT_H=$2; shift 2 ;;
        --samples) SAMPLES=$2; shift 2 ;;
        --resources) RESOURCES=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
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
command -v taskset >/dev/null || die "taskset not found (needed to enforce the thread limit)"
CPU_SLOTS=$(python3 scripts/g2f_tool.py cpu-slots --jobs "$JOBS" --threads "$THREADS" 2>&1) \
    || die "$CPU_SLOTS"
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
FREE_GB=$(( $(df --output=avail -B1 "$RESDIR" | tail -1) / 1000000000 ))
say "free disk on $(df --output=target "$RESDIR" | tail -1 | tr -d ' '): ${FREE_GB} GB; reserve ${RESERVE_GB} GB; prepared HiFi FASTQ retention: ${RETAIN_MODE}"
MEM_AVAIL_GB=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo)
ARC_GB=$(awk '$1=="size"{print int($3/1e9)}' /proc/spl/kstat/zfs/arcstats 2>/dev/null || echo 0)
say "memory: ${MEM_AVAIL_GB} GB available + ${ARC_GB} GB ZFS cache (reclaimable); budget ${JOBS} x ${MEM_GB} GB"
(( JOBS * MEM_GB <= MEM_AVAIL_GB + ARC_GB )) || say "WARNING: memory budget exceeds available memory"
say "load average: $(cut -d' ' -f1-3 /proc/loadavg); threads requested: ${JOBS} x ${THREADS} = $((JOBS * THREADS))"
say "CPU blocks, one per running genotype (enforced with taskset): $CPU_SLOTS"

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

# Disk projection per genotype (scripts/g2f/diskplan.py) and the plan for this launch.
STATE="$LOGDIR/state"; mkdir -p "$STATE"
python3 scripts/g2f_tool.py discover --json 2>/dev/null \
    | python3 scripts/g2f_tool.py disk-estimates --out "$STATE/estimates.json" --only $SAMPLES \
    || die "disk estimates failed"
python3 - "$STATE/estimates.json" "$FREE_GB" "$RESERVE_GB" "$JOBS" <<'EOF' | tee -a "$LOGDIR/launcher.log"
import json, sys
est = json.load(open(sys.argv[1])); free, reserve, jobs = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
G = 1e9
print(f"{'genotype':10s} {'peak GB':>8s} {'kept GB (FASTQ kept)':>21s} {'kept GB (FASTQ deleted)':>24s}")
for s, e in est.items():
    print(f"{s:10s} {e['peak']/G:8.0f} {e['persist_keep']/G:21.0f} {e['persist_temp']/G:24.0f}")
temp = sum(e["persist_temp"] for e in est.values())
keep = sum(e["persist_keep"] for e in est.values())
trans = sorted((e["peak"] - e["persist_temp"] for e in est.values()), reverse=True)
worst = temp + sum(trans[:jobs])
print(f"final footprint: {temp/G:.0f} GB (all FASTQs deleted after use) to {keep/G:.0f} GB (all kept)")
print(f"worst-case peak with FASTQs deleted after use: {worst/G:.0f} GB; free {free:.0f} GB, reserve {reserve:.0f} GB")
if worst / G > free - reserve:
    print("WARNING: even with every prepared FASTQ deleted after use, the worst case exceeds free space "
          "minus reserve; the disk gate will make genotypes wait for space")
EOF

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
           "mem_gb_per_genotype": $MEM_GB, "cpu_slots": "$CPU_SLOTS".split(), "resources": "$RESOURCES", "samples": "$SAMPLES".split(),
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
    say "  taskset -c <CPU block> ${PER_GENO[*]} --config resources=$RESOURCES only_samples=<LINE>"
    exit 0
fi

# ------------------------------------------------------------------ execution
say "creating the shared validation-identity stamp"
snakemake identity_stamp "${COMMON[@]}" --cores 1 --config "resources=$RESOURCES" \
    >> "$LOGDIR/launcher.log" 2>&1 || die "could not create the identity stamp"

printf 'genotype\texit\tretain_prepared_hifi\tstarted\tended\tlog\n' > "$LOGDIR/status.tsv"
export LOGDIR STATE RESDIR RESOURCES THREADS MEM_GB RESERVE_GB JOBS RETAIN_MODE MAX_WAIT_H CPU_SLOTS

# Disk monitor: free space and results size every 10 minutes.
( set +e +o pipefail
  printf 'time\tfree_GB\tresults_GB\trunning\n' > "$LOGDIR/disk.tsv"
  while :; do
      f=$(( $(df --output=avail -B1 "$RESDIR" | tail -1) / 1000000000 ))
      r=$(du -s -B1 "$RESDIR" 2>/dev/null | awk '{printf "%d", $1/1e9}')
      run=$(cd "$STATE" && ls *.running 2>/dev/null | sed 's/.running//' | tr '\n' ',')
      printf '%s\t%s\t%s\t%s\n' "$(date '+%F %T')" "$f" "$r" "$run" >> "$LOGDIR/disk.tsv"
      (( f < RESERVE_GB )) && echo "[$(date '+%F %T')] WARNING: free disk ${f} GB is below the reserve (${RESERVE_GB} GB); no new genotype will start" >> "$LOGDIR/launcher.log"
      sleep 600
  done ) &
MONITOR=$!
trap 'kill $MONITOR 2>/dev/null' EXIT

say "starting: ${JOBS} genotypes at a time, ${THREADS} threads each; logs in $LOGDIR"
say "drain gracefully (finish running genotypes, start no new ones): touch $STATE/STOP"
printf '%s\n' $SAMPLES | parallel -j "$JOBS" --halt never --joblog "$LOGDIR/joblog.tsv" \
    bash scripts/run_one_genotype.sh {} {%} || true

python3 scripts/g2f_tool.py raw-release >> "$LOGDIR/launcher.log" 2>&1 || true
FAILED=$(awk -F'\t' 'NR>1 && $2!=0{print $1}' "$LOGDIR/status.tsv" | tr '\n' ' ')
if [[ -n "$FAILED" ]]; then
    say "NOT COMPLETED: $FAILED (see status.tsv and their logs; the others are complete). Fix, then"
    say "rerun this launcher: completed work is not repeated. Project summary not built."
    exit 1
fi
say "all genotypes completed; building the project summary table"
snakemake all "${COMMON[@]}" --cores 4 --config "resources=$RESOURCES" >> "$LOGDIR/launcher.log" 2>&1 \
    || die "project summary failed"
say "DONE. Summary: $(ls -t "$RESDIR"/summary/assembly_summary.*.tsv | head -1)"
say "Raw copies eligible for manual deletion: $RESDIR/raw_release_report.tsv"
exit 0
