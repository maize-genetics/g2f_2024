#!/bin/bash
#SBATCH --job-name=te_restore_ids
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00

# Step 06e: put the real sequence names back (EDTA ran on short names s00001, ...) and collect the
# final files in runs/<unit>/results/. Submitted by te_submit.sh after te_repeatmasker.

set -euo pipefail
: "${TE:?}" "${NAME:?}" "${SCRIPTS:?}"
echo "Date: $(date --iso-8601=seconds)   NAME=$NAME"
python3 "$SCRIPTS/te_restore_ids.py" --te "$TE" "$NAME"
echo "$(date --iso-8601=seconds) done"
