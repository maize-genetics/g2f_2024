#!/usr/bin/env bash
# Assembly validation, one sequential run: environment check -> all planned
# stages in order -> a new report per sample. Exits with the run's status
# (non-zero if any stage did not validate); the report is written either way
# and shows which stages are validated.
#
#   bash scripts/run_asmqc.sh CONFIG [run options, e.g. --stages ... --samples ...]
#
# Start it detached so the session ends with the run, e.g.
#   run_detached asmqc LOGDIR -- bash scripts/run_asmqc.sh config/local/asmqc.yaml
# Python: needs numpy, pysam, matplotlib, networkx, PyYAML (ASMQC_PYTHON overrides).
set -euo pipefail
[[ $# -ge 1 ]] || { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
CONFIG=$1; shift
cd "$(dirname "$0")/.."
PY=${ASMQC_PYTHON:-python3}
echo "== $(date '+%F %T') environment check"
"$PY" scripts/g2f_asmqc.py -c "$CONFIG" check-env
echo "== $(date '+%F %T') plan"
"$PY" scripts/g2f_asmqc.py -c "$CONFIG" plan "$@"
echo "== $(date '+%F %T') run"
rc=0
"$PY" scripts/g2f_asmqc.py -c "$CONFIG" run "$@" || rc=$?
echo "== $(date '+%F %T') report"
SAMPLES=()
for ((i = 1; i <= $#; i++)); do
    if [[ ${!i} == --samples ]]; then j=$((i + 1)); SAMPLES=(--samples "${!j}"); fi
done
"$PY" scripts/g2f_asmqc.py -c "$CONFIG" report "${SAMPLES[@]}" || rc=${rc:-1}
echo "== $(date '+%F %T') finished, run exit status $rc"
exit "$rc"
