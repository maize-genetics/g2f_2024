#!/bin/bash
#SBATCH --job-name=te_lai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 06f: LTR Assembly Index (LAI) of the nuclear genome, one array task per line of
# PROJECT_DIR/units.tsv. Uses the LAI program of LTR_retriever in the EDTA 2.3.0 image.
#
# Reads (TE/runs/<unit>/results/, real sequence names):
#   <unit>.nuclear.softmasked.fa      nuclear genome sequence
#   <unit>.LTR.pass.list              intact LTR retrotransposons
#   <unit>.nuclear.RepeatMasker.out   all TE hits
# Writes:
#   TE/runs/<unit>/lai/<unit>.nuclear.RepeatMasker.out.LAI   (whole_genome row = LAI of the nuclear genome)
#
# LAI can take more than a day per genome; the limit here is 4 days.
#
# Usage (from the repository folder), N = number of lines in units.tsv:
#   source config.sh && N=$(wc -l < "$PROJECT_DIR/units.tsv")
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1-$N 06_te_annotation/te_lai.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
[[ -f ${EDTA230_SIF:-} ]] || { echo "Set EDTA230_SIF in config.sh"; exit 1; }
command -v apptainer >/dev/null || { echo "apptainer not found"; exit 1; }

TE=$PROJECT_DIR/te_annotation
U=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJECT_DIR/units.tsv" | cut -f1)
[[ -n $U ]] || { echo "No unit on line ${SLURM_ARRAY_TASK_ID} of $PROJECT_DIR/units.tsv"; exit 1; }
R=$TE/runs/$U/results
W=$TE/runs/$U/lai

echo "$(date --iso-8601=seconds)  LAI  $U  job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}  $(hostname)"
for f in "$R/$U.nuclear.softmasked.fa" "$R/$U.LTR.pass.list" "$R/$U.nuclear.RepeatMasker.out"; do
  [[ -s $f ]] || { echo "Failed. missing $f"; exit 1; }
done
echo "intact LTRs in pass list: $(grep -vc '^#' "$R/$U.LTR.pass.list")"

mkdir -p "$W"
cd "$W"
ln -sf "$R/$U.nuclear.softmasked.fa" genome.fa
ln -sf "$R/$U.LTR.pass.list" "$U.LTR.pass.list"
ln -sf "$R/$U.nuclear.RepeatMasker.out" "$U.nuclear.RepeatMasker.out"

apptainer exec --cleanenv --bind "$TE" "$EDTA230_SIF" \
  LAI -genome genome.fa -intact "$U.LTR.pass.list" -all "$U.nuclear.RepeatMasker.out" -t "${SLURM_CPUS_PER_TASK:-1}"

out=$(ls -t "$W"/*.LAI 2>/dev/null | head -1 || true)
[[ -n $out ]] || { echo "Failed. no .LAI file written"; exit 1; }
echo "LAI file: $out"
grep -E "^Chr|whole_genome" "$out" || true
echo "$(date --iso-8601=seconds) done"
