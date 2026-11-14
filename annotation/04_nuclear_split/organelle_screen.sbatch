#!/bin/bash
#SBATCH --job-name=organelle_screen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 04a: label every scaffold of each assembly as nuclear, plastid or mitochondrion, by aligning
# the scaffolds to the maize plastid and mitochondrial genomes with minimap2.
#
# Reads:  PROJECT_DIR/units.tsv, PROJECT_DIR/inputs/<unit>.scaffolds.fa, ORGANELLE_FASTA
#         PROJECT_DIR/scaffold_classes/review_overrides.tsv   (optional; your decisions for 'review' scaffolds)
# Writes: PROJECT_DIR/scaffold_classes/<unit>.vs_organelle.paf.gz
#         PROJECT_DIR/scaffold_classes/<unit>.scaffold_classification.tsv
#         PROJECT_DIR/scaffold_classes/scaffold_summary.tsv
#
# Usage (from the repository folder):
#   source config.sh && sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" 04_nuclear_split/organelle_screen.sbatch
# Rerun after adding decisions to review_overrides.tsv (the alignments are reused).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
: "${ORGANELLE_FASTA:?set ORGANELLE_FASTA in config.sh}"
MINIMAP2=${MINIMAP2:-minimap2}
HERE=${SLURM_SUBMIT_DIR:-.}/04_nuclear_split

BASE=$PROJECT_DIR
OUT=$BASE/scaffold_classes
OVR=$OUT/review_overrides.tsv
mkdir -p "$OUT"
[[ -f $OVR ]] || printf "unit\tscaffold\tfinal_class\tconfidence\treason\n" > "$OVR"

echo "minimap2 $($MINIMAP2 --version)"
SUMMARY=$OUT/scaffold_summary.tsv
printf "unit\tnuclear_scaffolds\tnuclear_Mb\tplastid_scaffolds\tplastid_Mb\tmito_scaffolds\tmito_Mb\tunresolved_n\tunresolved\n" > "$SUMMARY"

while IFS=$'\t' read -r UNIT _; do
  [[ -n $UNIT ]] || continue
  echo "== $UNIT  $(date)"
  FA=$BASE/inputs/${UNIT}.scaffolds.fa
  PAF=$OUT/${UNIT}.vs_organelle.paf.gz
  if [[ ! -s $PAF ]]; then
    $MINIMAP2 -x asm5 -c -t "${SLURM_CPUS_PER_TASK:-1}" "$ORGANELLE_FASTA" "$FA" 2> "$OUT/${UNIT}.minimap2.log" \
      | gzip > "$PAF.partial"
    mv "$PAF.partial" "$PAF"
  fi
  python3 "$HERE/classify_scaffolds.py" "$UNIT" "$FA" "$PAF" "$OUT/${UNIT}.scaffold_classification.tsv" "$OVR" >> "$SUMMARY"
done < "$BASE/units.tsv"

echo
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
echo
echo "Scaffolds still to decide (final_class = review; add a decision to $OVR and rerun):"
awk -F'\t' 'FNR>1 && $8=="review"{n=split(FILENAME,a,"/"); print a[n]"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6}' \
  "$OUT"/*.scaffold_classification.tsv || true
echo "Done $(date)"
