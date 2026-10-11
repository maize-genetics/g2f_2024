#!/usr/bin/env bash
# Run the lightweight test suite (synthetic data only; ~1 minute).
#   bash tests/run_tests.sh            # all tests
#   bash tests/run_tests.sh test_units # one module
# On BioHPC, Snakemake 9 and PyYAML come from the Python 3.12 module.
set -euo pipefail
cd "$(dirname "$0")"
if command -v module >/dev/null 2>&1; then module load python/3.12.7 >/dev/null 2>&1 || true; fi
export TMPDIR="${TMPDIR:-/workdir/$USER/tmp}"
mkdir -p "$TMPDIR"
if [[ $# -eq 0 ]]; then set -- test_discovery test_units test_workflow; fi
python3 -m unittest -v "$@"
