#!/bin/bash
#SBATCH --job-name=helixer
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 01b: Helixer gene prediction on each whole genome assembly, one array task per line of
# PROJECT_DIR/units.tsv (written by helixer_prep_inputs.sbatch).
# Helixer 0.3.3, lineage land_plant (default model land_plant_v0.3_a_0080, default window size).
#
# Writes: PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.raw.gff3
#         PROJECT_DIR/raw/<unit>/<unit>.run_info.txt       versions, model md5, parameters
#         PROJECT_DIR/raw/<unit>/<unit>.summary_gff3.txt   feature counts
#
# Usage (from the repository folder), N = number of lines in units.tsv:
#   source config.sh && N=$(wc -l < "$PROJECT_DIR/units.tsv")
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 01_helixer/helixer_array.sbatch          # test one
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N%4 01_helixer/helixer_array.sbatch     # the rest
# Add  --mail-type=FAIL,END --mail-user="$EMAIL"  for notifications.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
HELIXER=${HELIXER:-Helixer.py}          # how to call Helixer on your system (see config.sh)

BASE=$PROJECT_DIR
UNITS=$BASE/units.tsv
CHECK=$BASE/inputs/prep_check.tsv

UNIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$UNITS" | cut -f1)
[[ -n $UNIT ]] || { echo "No unit on line ${SLURM_ARRAY_TASK_ID} of $UNITS"; exit 1; }
FASTA=$BASE/inputs/${UNIT}.scaffolds.fa

# only run on inputs that passed the prep check
if ! awk -F'\t' -v u="$UNIT" '$1==u && $3=="OK"{f=1} END{exit !f}' "$CHECK"; then
  echo "$UNIT did not pass the prep check ($CHECK); not running."; exit 1
fi

OUT=$BASE/raw/$UNIT
GFF=$OUT/${UNIT}.scaffolds.helixer.raw.gff3
mkdir -p "$OUT"
# don't overwrite a finished run if the task is resubmitted
[[ -s $GFF ]] && { echo "$GFF already exists; skipping $UNIT."; exit 0; }

# record what was run, for the methods
INFO=$OUT/${UNIT}.run_info.txt
{
  echo "unit: $UNIT"
  echo "input: $FASTA"
  echo "job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $(hostname)"
  echo "start: $(date)"
  echo "gpu: $(nvidia-smi -L 2>/dev/null || echo none found)"
  echo "helixer: $($HELIXER --version 2>/dev/null | tail -n 1)"   # expect Helixer.py 0.3.3
  echo "params: --lineage land_plant (default window size and model), --species $UNIT"
} | tee "$INFO"

# write to a temporary name so a crash never leaves a GFF3 that looks finished
$HELIXER \
  --lineage land_plant \
  --fasta-path "$FASTA" \
  --species "$UNIT" \
  --gff-output-path "$GFF.partial"
mv "$GFF.partial" "$GFF"

echo "end: $(date)" | tee -a "$INFO"

# quick summary
{
  for t in gene mRNA exon CDS five_prime_UTR three_prime_UTR; do
    printf "%s\t%s\n" "$t" "$(awk -v t="$t" '$3==t' "$GFF" | wc -l)"
  done
  printf "sequences_with_genes\t%s\n" "$(awk '$3=="gene"{print $1}' "$GFF" | sort -u | wc -l)"
} | tee "$OUT/${UNIT}.summary_gff3.txt"

echo "Helixer $UNIT complete"
