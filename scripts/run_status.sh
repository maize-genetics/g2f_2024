#!/usr/bin/env bash
# Read-only status of the full-scale run: which launcher run is current, the
# state of every genotype, the running genotype's progress, disk, CPU, memory.
# Usage (from anywhere):  bash scripts/run_status.sh [logs/assembly_<stamp>]
# States: done (completed, exit 0) | running | failed | interrupted (launcher
# gone: crash/kill) | not_started (queue ended first) | pending (not reached yet).
set -uo pipefail
cd "$(dirname "$0")/.."
LOGROOT=$(realpath ../logs)
RUN=${1:-$(ls -d "$LOGROOT"/assembly_2* 2>/dev/null | grep -v '\.dryrun$' | sort | tail -1)}
[[ -d $RUN ]] || { echo "no launcher run found under $LOGROOT"; exit 1; }
if [[ -n ${PIPELINE_CONFIG:-} ]]; then TOPCFG=$PIPELINE_CONFIG
elif [[ -f config/local/config.yaml ]]; then TOPCFG=config/local/config.yaml
else TOPCFG=config/config.yaml; fi
RESDIR=$(python3 -c "import os,yaml;print(os.path.normpath(os.path.join('$(pwd)', yaml.safe_load(open('$TOPCFG'))['results_dir'])))")

if flock -n "$LOGROOT/launcher.lock" true 2>/dev/null; then LIVE="no launcher running"
else LIVE="launcher RUNNING"; fi
echo "run: $RUN ($LIVE)"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print("  started", d["started"], "| jobs", d["jobs"], "x", d["threads"], "threads | cpus", " ".join(d["cpu_slots"]), "| on failure", d.get("on_failure", "continue"))' "$RUN/launch.json" 2>/dev/null
[[ -e $RUN/state/STOP ]] && echo "  STOP file present: no further genotype will start"

echo; printf '%-10s %-12s %s\n' genotype state detail
for s in $(python3 -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["samples"]))' "$RUN/launch.json"); do
    st=pending; detail=""
    for k in done failed interrupted running; do [[ -e $RUN/state/$s.$k ]] && st=$k; done
    row=$(awk -F'\t' -v s="$s" 'NR>1 && $1==s' "$RUN/status.tsv" 2>/dev/null | tail -1)
    if [[ -n $row ]]; then
        IFS=$'\t' read -r _ ex _ t0 t1 _ <<< "$row"
        [[ $st == pending && $ex == not_started ]] && st=not_started
        detail="exit=$ex $t0 -> $t1"
    fi
    if [[ $st == running && $LIVE == "no launcher running" ]]; then
        st=interrupted; detail="launcher gone (crash/kill); marked on the next launch"
    elif [[ $st == running ]]; then
        detail="$(grep -E '^[0-9]+ of [0-9]+ steps' "$RUN/$s.log" 2>/dev/null | tail -1)"
        h=$(ls -t "$RESDIR/$s"/assembly/*/"$s".asm.stderr.log 2>/dev/null | head -1)
        if [[ -n $h ]] && pgrep -u "$(id -u)" -f "hifiasm -o $s.asm" >/dev/null; then
            detail="$detail | hifiasm: $(grep -E 'corrected reads for round|bin files have been written|Index has been written|Writing .* GFA|Real time' "$h" | tail -1 | cut -c1-80)"
        fi
    fi
    printf '%-10s %-12s %s\n' "$s" "$st" "$detail"
done

echo
echo "disk: $(df -h --output=avail,pcent "$RESDIR" | tail -1 | awk '{print $1 " free (" $2 " used)"}') on $(df --output=target "$RESDIR" | tail -1); last entries of disk.tsv:"
tail -3 "$RUN/disk.tsv" 2>/dev/null | sed 's/^/  /'
echo "load: $(cut -d' ' -f1-3 /proc/loadavg) on $(nproc --all) CPUs; memory: $(free -g | awk '/Mem:/{print $7 " GB available of " $2}')"
echo "our processes (snakemake/hifiasm/busco/merqury/meryl):"
ps -u "$(id -u)" -o pid,etime,pcpu,rss,args --sort=-pcpu 2>/dev/null \
    | awk 'NR>1 && /snakemake|hifiasm|busco|merqury|meryl|pbtk|bam2fastq/ && !/awk/ {printf "  %s %s %s%% %dGB %s\n", $1, $2, $3, $4/1048576, substr($0, index($0,$5), 90)}' | head -8
echo; echo "launcher log (last 5 lines):"; tail -5 "$RUN/launcher.log" | sed 's/^/  /'
