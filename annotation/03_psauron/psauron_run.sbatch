#!/bin/bash
#SBATCH --job-name=psauron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 03a: PSAURON 1.1.3 score for every Helixer CDS of the whole genome, on CPU (-c).
# One array task per line of PROJECT_DIR/units.tsv.
#
# Reads:  PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.cds.fa      (from step 02)
# Writes: PROJECT_DIR/psauron/<unit>.helixer.cds.psauron.csv           one row per transcript
#         (the overall PSAURON score of the genome is printed in the log and in the CSV's first lines)
#
# Usage (from the repository folder), N = number of lines in units.tsv:
#   source config.sh && N=$(wc -l < "$PROJECT_DIR/units.tsv")
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 03_psauron/psauron_run.sbatch       # test one
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N 03_psauron/psauron_run.sbatch    # the rest

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
PSAURON=${PSAURON:-psauron}

BASE=$PROJECT_DIR
UNIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE/units.tsv" | cut -f1)
[[ -n $UNIT ]] || { echo "No unit on line ${SLURM_ARRAY_TASK_ID} of $BASE/units.tsv"; exit 1; }
CDS=$BASE/raw/$UNIT/${UNIT}.scaffolds.helixer.cds.fa
[[ -s $CDS ]] || { echo "missing $CDS; run 02_busco/busco_extract_seqs.sbatch first"; exit 1; }

mkdir -p "$BASE/psauron"
OUT=$BASE/psauron/${UNIT}.helixer.cds.psauron.csv
[[ -s $OUT ]] && { echo "$OUT already exists; skipping $UNIT."; exit 0; }

echo "unit: $UNIT   input: $CDS ($(grep -c '^>' "$CDS") sequences)   $(date)   $(hostname)"

# write to a temporary name so a crash never leaves a CSV that looks finished
$PSAURON -c -i "$CDS" -o "$OUT.part"
mv "$OUT.part" "$OUT"

grep -m1 -i "psauron score" "$OUT" || true
echo "done $(date): $OUT"
