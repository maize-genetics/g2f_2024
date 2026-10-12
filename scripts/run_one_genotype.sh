#!/usr/bin/env bash
# One genotype of the full-scale run (called by scripts/run_all_assemblies.sh
# through GNU parallel; not meant to be run by hand).
#   1. wait until the disk gate allows the start (scripts/g2f/diskplan.py);
#      the gate also decides whether the prepared HiFi FASTQ is retained
#   2. run the workflow for this genotype (target per_sample)
#   3. record the outcome and refresh the raw-copy release report (no deletion)
# Arguments: genotype, GNU parallel job slot (1..JOBS).
# Environment (exported by the launcher): LOGDIR STATE RESDIR RESOURCES THREADS
# MEM_GB RESERVE_GB JOBS RETAIN_MODE MAX_WAIT_H CPU_SLOTS (one CPU list per
# slot; the genotype's Snakemake and every tool it runs are pinned to it)
set -uo pipefail
S=$1
SLOT=${2:-}
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> "$LOGDIR/launcher.log"; }
status() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$S" "$1" "$2" "${START:-}" "$(date '+%F %T')" \
    "$LOGDIR/$S.log" >> "$LOGDIR/status.tsv"; }

PIN=(); CPUS=unpinned
if [[ -n "${CPU_SLOTS:-}" && -n "$SLOT" ]]; then
    read -ra blocks <<< "$CPU_SLOTS"
    CPUS=${blocks[$((SLOT - 1))]:-}
    if [[ -z "$CPUS" ]]; then
        log "$S: no CPU block for slot $SLOT (CPU_SLOTS=$CPU_SLOTS)"; status cpu_slot_error -; exit 6
    fi
    PIN=(taskset -c "$CPUS")
fi

waited=0
while :; do
    if [[ -e "$STATE/STOP" ]]; then
        log "$S not started: STOP file present ($STATE/STOP)"; status stopped -; exit 3
    fi
    out=$(python3 scripts/g2f_tool.py disk-gate --sample "$S" --state-dir "$STATE" \
          --results "$RESDIR" --reserve-gb "$RESERVE_GB" --jobs "$JOBS" --mode "$RETAIN_MODE")
    rc=$?
    echo "$out" >> "$LOGDIR/gate.jsonl"
    [[ $rc -eq 0 ]] && break
    if [[ $rc -ne 10 && $rc -ne 11 ]]; then
        log "$S: disk gate error (exit $rc): $out"; status gate_error -; exit 4
    fi
    if (( waited >= MAX_WAIT_H * 3600 )); then
        log "$S: not started after ${MAX_WAIT_H} h waiting for disk space"; status disk_wait_timeout -; exit 5
    fi
    log "$S waiting for disk space ($( [[ $rc -eq 11 ]] && echo 'does not fit even alone; ' )$out)"
    sleep 600; waited=$((waited + 600))
done
RETAIN=$(python3 -c 'import json,sys; print(str(json.loads(sys.argv[1])["retain_prepared_hifi"]).lower())' "$out")
START=$(date '+%F %T')
log "START $S (retain_prepared_hifi=$RETAIN, cpus=$CPUS): $out"

"${PIN[@]}" snakemake per_sample --profile profiles/biohpc --rerun-triggers mtime --cores "$THREADS" \
    --resources "mem_mb=$((MEM_GB * 1000))" --keep-going --rerun-incomplete \
    --config "resources=$RESOURCES" "only_samples=$S" "retain_prepared_hifi=$RETAIN" \
    > "$LOGDIR/$S.log" 2>&1
rc=$?
mv "$STATE/$S.running" "$STATE/$S.$([[ $rc -eq 0 ]] && echo done || echo failed)"
python3 scripts/g2f_tool.py raw-release --only "$S" >> "$LOGDIR/$S.raw_release.txt" 2>&1 || true
du -sh "$RESDIR/$S" 2>/dev/null | sed "s|^|[$(date '+%F %T')] results size: |" >> "$LOGDIR/launcher.log"
log "END $S exit=$rc"
status "$rc" "$RETAIN"
exit $rc
