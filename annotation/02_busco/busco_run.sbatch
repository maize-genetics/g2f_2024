#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 02b: BUSCO 6.1.0 in protein mode with poales_odb12, on the proteins of every gene in the
# whole genome. One array task per line of PROJECT_DIR/units.tsv. Runs offline from a lineage
# folder downloaded beforehand (BUSCO_LINEAGE in config.sh).
#
# Reads:  PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.proteins.fa   (from busco_extract_seqs.sbatch)
# Writes: PROJECT_DIR/busco/runs/<unit>/                                BUSCO's own folder, incl. short_summary*.txt
#
# Usage (from the repository folder), N = number of lines in units.tsv:
#   source config.sh && N=$(wc -l < "$PROJECT_DIR/units.tsv")
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 02_busco/busco_run.sbatch       # test one
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N 02_busco/busco_run.sbatch    # the rest

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
: "${BUSCO_LINEAGE:?set BUSCO_LINEAGE in config.sh (downloaded poales_odb12 folder)}"
BUSCO=${BUSCO:-busco}

BASE=$PROJECT_DIR
RUNS=$BASE/busco/runs
UNIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE/units.tsv" | cut -f1)
[[ -n $UNIT ]] || { echo "No unit on line ${SLURM_ARRAY_TASK_ID} of $BASE/units.tsv"; exit 1; }
PROT=$BASE/raw/$UNIT/${UNIT}.scaffolds.helixer.proteins.fa
[[ -s $PROT ]] || { echo "missing $PROT; run busco_extract_seqs.sbatch first"; exit 1; }

mkdir -p "$RUNS"
# don't rerun a unit that already has a result
if ls "$RUNS/$UNIT"/short_summary*.txt >/dev/null 2>&1; then
  echo "$UNIT already has a BUSCO result; skipping."; exit 0
fi

echo "unit: $UNIT   input: $PROT   $(date)   $(hostname)"
$BUSCO --version
grep -E "^(name|dataset_version|number_of_BUSCOs)=" "$BUSCO_LINEAGE/dataset.cfg" || true

$BUSCO \
  -i "$PROT" \
  -m proteins \
  -l "$BUSCO_LINEAGE" \
  --offline \
  -c "${SLURM_CPUS_PER_TASK:-1}" \
  -o "$UNIT" \
  --out_path "$RUNS"

echo "== $UNIT result"
grep -E "C:[0-9.]+%" "$RUNS/$UNIT"/short_summary*.txt
echo "done $(date)"
