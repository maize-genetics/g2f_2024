#!/usr/bin/env bash
# Full-scale launcher: assemble and QC every genotype with the existing
# Snakemake workflow, ONE genotype at a time by default; --jobs N runs N at a
# time. Defaults come from the `launcher:` section of the top-level
# configuration ($PIPELINE_CONFIG, else config/local/config.yaml, else
# config/config.yaml); command-line options override them.
#
# Each genotype runs as its own Snakemake process (target `per_sample`,
# `--config only_samples=<LINE>`), so logs are per genotype and only one
# genotype holds large intermediates at a time. The next genotype starts only
# after the previous one completed (exit 0); by default the first failure,
# STOP or refusal ends the queue (--on-failure continue: carry on). Rerunning
# the launcher resumes: completed stages are not repeated (Snakemake),
# completed genotypes finish immediately, and genotypes left `running` by a
# crashed launcher are marked `interrupted` in their old log directory.
#
# Never more than N genotypes run: GNU parallel -j N, a hard limit in the disk
# gate (scripts/g2f/diskplan.py), and a project-wide lock (logs/launcher.lock)
# that refuses a second launcher, or a launch while Snakemake/hifiasm processes
# are working in this pipeline or results directory. Each genotype's Snakemake,
# and every tool it runs, is pinned (taskset) to its own block of T CPUs, one
# per physical core.
#
# Run from pipeline/. Review first:
#   bash scripts/run_all_assemblies.sh --dry-run
# Then start it in a detached session that closes when the launcher ends
# (docs/reference.md, Running all samples):
#   S=$(date +%Y%m%d_%H%M%S)
#   screen -dmS g2f_run bash -c "set -o pipefail; bash scripts/run_all_assemblies.sh 2>&1 | tee ../logs/launcher.$S.out; echo \$? > ../logs/launcher.$S.exit"
# Monitor:  bash scripts/run_status.sh
#
# Options:
#   --jobs N            genotypes at a time (default 1)
#   --threads T         threads per genotype (default 32). Enforced: each running
#                       genotype is pinned (taskset) to its own block of T CPUs,
#                       so the run never uses more than jobs x threads CPUs
#   --mem-gb G          memory budget per genotype for Snakemake scheduling (default 128)
#   --reserve-gb R      free disk always kept in reserve by the disk gate (default 100)
#   --retain-hifi M     prepared HiFi FASTQ retention: keep (default; nothing is
#                       deleted), auto (the disk gate may make it temporary per
#                       genotype), or temp (always delete after use). With auto or
#                       temp, keep the raw BAMs: they become the only local copy.
#   --max-wait-h H      a genotype waiting this long for disk space fails (default 1)
#   --on-failure M      stop (default): a genotype that does not complete ends the
#                       queue; continue: start the next genotype anyway
#   --samples "A B"     only these genotypes, in this order (default: all discovered;
#                       assembled first, then hifiasm checkpoints, then smallest HiFi)
#   --resources FILE    resources file (default config/resources/default.yaml)
#   --profile DIR       Snakemake profile (default profiles/local)
#   --dry-run           preflight + `snakemake -n` per genotype; executes nothing
# (The defaults above are those of the tracked config/config.yaml.)
set -euo pipefail

cd "$(dirname "$0")/.."
PIPE=$(pwd)
# Environment modules (e.g. Cornell BioHPC); a no-op where the module does not exist.
if command -v module >/dev/null 2>&1; then module load python/3.12.7 >/dev/null 2>&1 || true; fi
JOBS=1; THREADS=32; MEM_GB=128; RESERVE_GB=100; RETAIN_MODE=keep; MAX_WAIT_H=1
SAMPLES=""; DRY=0; ON_FAILURE=stop
RESOURCES=config/resources/default.yaml; PROFILE=profiles/local
if [[ -n ${PIPELINE_CONFIG:-} ]]; then TOPCFG=$PIPELINE_CONFIG
elif [[ -f config/local/config.yaml ]]; then TOPCFG=config/local/config.yaml
else TOPCFG=config/config.yaml; fi
if [[ -f $TOPCFG ]]; then
    DEFAULTS=$(python3 - "$TOPCFG" <<'PYEOF'
import shlex, sys, yaml
names = {"jobs": "JOBS", "threads": "THREADS", "mem_gb": "MEM_GB", "reserve_gb": "RESERVE_GB",
         "retain_hifi": "RETAIN_MODE", "max_wait_h": "MAX_WAIT_H", "on_failure": "ON_FAILURE",
         "resources": "RESOURCES", "profile": "PROFILE"}
with open(sys.argv[1]) as fh:
    d = (yaml.safe_load(fh) or {}).get("launcher") or {}
bad = sorted(set(d) - set(names))
if bad:
    sys.exit(f"{sys.argv[1]}: unknown launcher setting(s): {', '.join(bad)}")
for k, v in d.items():
    print(f"{names[k]}={shlex.quote(str(v))}")
PYEOF
) || exit 2
    eval "$DEFAULTS"
fi
while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs) JOBS=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --mem-gb) MEM_GB=$2; shift 2 ;;
        --reserve-gb) RESERVE_GB=$2; shift 2 ;;
        --retain-hifi) RETAIN_MODE=$2; shift 2 ;;
        --max-wait-h) MAX_WAIT_H=$2; shift 2 ;;
        --on-failure) ON_FAILURE=$2; shift 2 ;;
        --samples) SAMPLES=$2; shift 2 ;;
        --resources) RESOURCES=$2; shift 2 ;;
        --profile) PROFILE=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) sed -n '2,52p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

PIPELINE_CONFIG=$(realpath -m "$TOPCFG"); export PIPELINE_CONFIG
unset OMP_NUM_THREADS   # rules that need it set it themselves
STAMP=$(date +%Y%m%d_%H%M%S)
LOGDIR=$(realpath -m "$PIPE/../logs/assembly_$STAMP")
[[ $DRY -eq 1 ]] && LOGDIR="$LOGDIR.dryrun"
mkdir -p "$LOGDIR"
say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOGDIR/launcher.log"; }
die() { say "ERROR: $*"; exit 1; }

[[ $ON_FAILURE == stop || $ON_FAILURE == continue ]] || die "--on-failure must be stop or continue"

# ------------------------------------------------------------------ exclusivity
# One launcher per project: the lock is held until this process exits and is
# released by the kernel however it ends (crash, kill, reboot), so it never
# goes stale. Dry runs execute nothing and do not take it.
LOGROOT=$(dirname "$LOGDIR")
if [[ $DRY -eq 0 ]]; then
    exec 9>"$LOGROOT/launcher.lock"
    flock -n 9 || die "another launcher is running for this project (lock $LOGROOT/launcher.lock)"
    # Also refuse while a Snakemake or hifiasm process of this user works in
    # this pipeline or results directory (e.g. one started by hand, or by an
    # older launcher without the lock). Other projects' runs do not count.
    PAT='^\S+ \S*/snakemake |^\S*snakemake |^\S*hifiasm '
    RES_EARLY=$(python3 -c 'import os,sys,yaml;print(os.path.normpath(os.path.join(sys.argv[1], (yaml.safe_load(open(sys.argv[2])) or {}).get("results_dir", "../results"))))' "$PIPE" "$PIPELINE_CONFIG" 2>/dev/null \
        || echo "$PIPE/../results")
    RES_EARLY=$(realpath -m "$RES_EARLY")
    BUSY=""
    for pid in $(pgrep -u "$(id -u)" -f "$PAT" || true); do
        cwd=$(readlink "/proc/$pid/cwd" 2>/dev/null) || continue
        case "$cwd/" in "$PIPE"/*|"$RES_EARLY"/*) BUSY+="$pid " ;; esac
    done
    if [[ -n $BUSY ]]; then
        ps -o pid=,args= -p "$(echo $BUSY | tr ' ' ,)" | tee -a "$LOGDIR/launcher.log" >&2
        die "Snakemake/hifiasm processes are still running in this project (listed above); stop them or wait"
    fi
    # Genotypes left `running` by an earlier launcher can only be interrupted
    # runs now (we hold the lock): mark them, so status is never ambiguous.
    for f in "$LOGROOT"/assembly_*/state/*.running; do
        [[ -e $f && $f != "$LOGDIR"/* ]] || continue
        g=$(basename "$f" .running); old=$(dirname "$(dirname "$f")")
        mv "$f" "${f%.running}.interrupted"
        printf '%s\tinterrupted\t-\t-\t%s\t%s\n' "$g" "$(date '+%F %T')" "$old/$g.log" >> "$old/status.tsv"
        say "marked $g in $(basename "$old") as interrupted (its launcher is gone)"
    done
fi

# ------------------------------------------------------------------ preflight
say "configuration: $PIPELINE_CONFIG; profile: $PROFILE"
say "pipeline: $PIPE (git $(git rev-parse --short HEAD 2>/dev/null || echo '?')$(git diff --quiet 2>/dev/null || echo ', UNCOMMITTED CHANGES'))"
[[ -f "$RESOURCES" ]] || die "resources file not found: $RESOURCES"
[[ -f "$PROFILE/config.yaml" ]] || die "Snakemake profile not found: $PROFILE"
command -v snakemake >/dev/null || die "snakemake (>= 8) not found on PATH"
SMK_VER=$(snakemake --version); [[ ${SMK_VER%%.*} -ge 8 ]] || die "Snakemake >= 8 required (found $SMK_VER)"
command -v parallel >/dev/null || die "GNU parallel not found"
command -v taskset >/dev/null || die "taskset not found (needed to enforce the thread limit)"
CPU_SLOTS=$(python3 scripts/g2f_tool.py cpu-slots --jobs "$JOBS" --threads "$THREADS" 2>&1) \
    || die "$CPU_SLOTS"
python3 - "$PIPE" <<'EOF' || die "tool check failed (see above)"
import os, shlex, sys, yaml
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
from g2f.config import load_config
with open(os.environ["PIPELINE_CONFIG"]) as fh:
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

RESDIR=$(python3 -c "import os,yaml;print(os.path.normpath(os.path.join('$PIPE', yaml.safe_load(open('$PIPELINE_CONFIG'))['results_dir'])))")
mkdir -p "$RESDIR"
FREE_GB=$(( $(df --output=avail -B1 "$RESDIR" | tail -1) / 1000000000 ))
say "free disk on $(df --output=target "$RESDIR" | tail -1 | tr -d ' '): ${FREE_GB} GB; reserve ${RESERVE_GB} GB; prepared HiFi FASTQ retention: ${RETAIN_MODE}"
MEM_AVAIL_GB=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo)
ARC_GB=$(awk '$1=="size"{print int($3/1e9)}' /proc/spl/kstat/zfs/arcstats 2>/dev/null || echo 0)
say "memory: ${MEM_AVAIL_GB} GB available + ${ARC_GB} GB ZFS cache (reclaimable); budget ${JOBS} x ${MEM_GB} GB"
(( JOBS * MEM_GB <= MEM_AVAIL_GB + ARC_GB )) || say "WARNING: memory budget exceeds available memory"
say "load average: $(cut -d' ' -f1-3 /proc/loadavg); threads requested: ${JOBS} x ${THREADS} = $((JOBS * THREADS))"
say "CPU blocks, one per running genotype (enforced with taskset): $CPU_SLOTS"

# Disk projection per genotype (scripts/g2f/diskplan.py) and the plan for this launch.
STATE="$LOGDIR/state"; mkdir -p "$STATE"
python3 scripts/g2f_tool.py discover --config "$PIPELINE_CONFIG" --json 2>/dev/null \
    | python3 scripts/g2f_tool.py disk-estimates --out "$STATE/estimates.json" ${SAMPLES:+--only $SAMPLES} \
    || die "disk estimates failed"
# Default order (diskplan.start_order): genotypes already assembled (QC only),
# then those with a hifiasm checkpoint, then the rest, smallest HiFi first, so
# assemblies finish as early as possible and the least work is exposed to an outage.
if [[ -z "$SAMPLES" ]]; then
    SAMPLES=$(python3 -c '
import json, os, sys
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
from g2f import diskplan
print(" ".join(diskplan.start_order(json.load(open(sys.argv[2])), sys.argv[3])))
' "$PIPE" "$STATE/estimates.json" "$RESDIR") || die "start order failed"
fi
[[ -n "$SAMPLES" ]] || die "no genotypes discovered"
say "genotypes (in start order): $SAMPLES"
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

COMMON=(--profile "$PROFILE" --rerun-triggers mtime)
say "whole-project dry run (discovery, configuration, DAG) ..."
# --rerun-incomplete as in the per-genotype runs: after an interruption the
# killed jobs are rerun, not reported as an error (resume must work)
snakemake per_sample "${COMMON[@]}" --cores "$THREADS" -n --rerun-incomplete --config "resources=$RESOURCES" \
    > "$LOGDIR/dryrun_all.log" 2>&1 \
    || die "dry run failed; see $LOGDIR/dryrun_all.log"
grep -A40 '^Job stats' "$LOGDIR/dryrun_all.log" | tee -a "$LOGDIR/launcher.log" || true

python3 - "$LOGDIR/launch.json" <<EOF
import json, subprocess, sys
git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True).stdout.strip()
json.dump({"started": "$STAMP", "dry_run": bool($DRY), "git_commit": git("rev-parse", "HEAD"),
           "git_dirty": bool(git("status", "--porcelain")), "jobs": $JOBS, "threads": $THREADS,
           "mem_gb_per_genotype": $MEM_GB, "on_failure": "$ON_FAILURE", "cpu_slots": "$CPU_SLOTS".split(), "resources": "$RESOURCES", "config": "$PIPELINE_CONFIG", "profile": "$PROFILE", "samples": "$SAMPLES".split(),
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
# A crash leaves Snakemake's directory lock behind; no Snakemake process of
# this user is alive (checked above, under the launcher lock), so it is stale.
say "removing any stale Snakemake lock left by an interrupted run"
snakemake --unlock "${COMMON[@]}" --cores 1 --config "resources=$RESOURCES" >> "$LOGDIR/launcher.log" 2>&1 \
    || die "snakemake --unlock failed"
say "creating the shared validation-identity stamp"
snakemake identity_stamp "${COMMON[@]}" --cores 1 --config "resources=$RESOURCES" \
    >> "$LOGDIR/launcher.log" 2>&1 || die "could not create the identity stamp"

printf 'genotype\texit\tretain_prepared_hifi\tstarted\tended\tlog\n' > "$LOGDIR/status.tsv"
export LOGDIR STATE RESDIR RESOURCES PROFILE THREADS MEM_GB RESERVE_GB JOBS RETAIN_MODE MAX_WAIT_H CPU_SLOTS

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

say "starting: ${JOBS} genotype(s) at a time, ${THREADS} threads each, on failure: ${ON_FAILURE}; logs in $LOGDIR"
say "drain gracefully (finish running genotypes, start no new ones): touch $STATE/STOP"
# stop: after the first genotype that does not complete, start no further one
# (running ones finish); continue: start every genotype regardless.
HALT=never; [[ $ON_FAILURE == stop ]] && HALT=soon,fail=1
printf '%s\n' $SAMPLES | parallel -j "$JOBS" --halt "$HALT" --joblog "$LOGDIR/joblog.tsv" \
    bash scripts/run_one_genotype.sh {} {%} || true
for s in $SAMPLES; do      # genotypes the queue never reached: pending, not failed
    awk -F'\t' -v s="$s" 'NR>1 && $1==s{f=1} END{exit !f}' "$LOGDIR/status.tsv" \
        || printf '%s\tnot_started\t-\t-\t-\t-\n' "$s" >> "$LOGDIR/status.tsv"
done

python3 scripts/g2f_tool.py raw-release --config "$PIPELINE_CONFIG" >> "$LOGDIR/launcher.log" 2>&1 || true
FAILED=$(awk -F'\t' 'NR>1 && $2!=0{print $1 "(" $2 ")"}' "$LOGDIR/status.tsv" | tr '\n' ' ')
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
