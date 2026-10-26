#!/usr/bin/env bash
# Run the lightweight test suite (synthetic data only; ~1 minute).
#   bash tests/run_tests.sh            # all tests
#   bash tests/run_tests.sh test_units # one module
# On BioHPC, Snakemake 9 and PyYAML come from the Python 3.12 module.
set -euo pipefail
cd "$(dirname "$0")"
# Assembly-validation tests need numpy/pysam/matplotlib (BioHPC: the default python3, before the module).
ASMQC_PY=${ASMQC_PYTHON:-$(command -v python3)}
if command -v module >/dev/null 2>&1; then module load python/3.12.7 >/dev/null 2>&1 || true; fi
export TMPDIR="${TMPDIR:-$(realpath -m ../../tmp)}"   # the project tmp/ (not the shared /tmp)
mkdir -p "$TMPDIR"
if [[ $# -eq 0 ]]; then set -- test_discovery test_units test_workflow test_fingerprint test_realtools test_release_fasta test_provenance_report test_workflow_stages; fi
if [[ $# -eq 1 && $1 == test_asmqc ]]; then exec "$ASMQC_PY" -m unittest -v test_asmqc; fi
# AnchorWave tests: the container tests run only when ANCHORWAVE_IMAGE names the .sif image.
if [[ $# -eq 1 && $1 == test_anchorwave ]]; then exec "$ASMQC_PY" -m unittest -v test_anchorwave; fi
python3 -m unittest -v "$@"
if [[ " $* " == *" test_units "* && " $* " == *" test_realtools "* ]]; then  # the full default suite
    "$ASMQC_PY" -m unittest -v test_asmqc
    "$ASMQC_PY" -m unittest -v test_anchorwave
fi
