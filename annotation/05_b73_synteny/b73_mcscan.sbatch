#!/bin/bash
#SBATCH --job-name=b73_mcscan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Step 05b: one-to-one syntenic orthologs between B73 (Zm00001eb.1 canonical) and one assembly,
# with MCScan (jcvi.compara.catalog ortholog): LAST on CDS -> C-score 0.99 filter (reciprocal best
# matches) -> synteny blocks -> liftover -> quota 1:1 screen. Then writes the ortholog tables
# (b73_table.py). One array task per line of PROJECT_DIR/units.tsv. Run b73_prep.py once first.
#
# Reads:  PROJECT_DIR/b73_synteny/runs/<unit>/B73.bed, B73.cds, <unit>.bed, <unit>.cds
# Writes: PROJECT_DIR/b73_synteny/runs/<unit>/   MCScan files (B73.<unit>.lifted.1x1.anchors, ...)
#                                                and <unit>.B73_orthologs*.tsv
#
# Usage (from the repository folder, with the jcvi environment active), N = lines in units.tsv:
#   source config.sh && N=$(wc -l < "$PROJECT_DIR/units.tsv")
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 05_b73_synteny/b73_mcscan.sbatch       # test one
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N 05_b73_synteny/b73_mcscan.sbatch    # the rest

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
JCVI_PYTHON=${JCVI_PYTHON:-python}
HERE=$(cd "${SLURM_SUBMIT_DIR:-.}/05_b73_synteny" && pwd)

BO=$PROJECT_DIR/b73_synteny
UNIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJECT_DIR/units.tsv" | cut -f1)
[[ -n $UNIT ]] || { echo "No unit on line ${SLURM_ARRAY_TASK_ID} of $PROJECT_DIR/units.tsv"; exit 1; }
RUN=$BO/runs/$UNIT
[[ -s $RUN/$UNIT.bed && -s $RUN/B73.bed ]] || { echo "inputs missing in $RUN; run b73_prep.py first"; exit 1; }

echo "unit: $UNIT   $(date)   $(hostname)"
$JCVI_PYTHON -c "import jcvi; print('jcvi', jcvi.__version__)"
lastal --version
# the 1:1 (quota) step needs an optimisation solver; stop now rather than produce a non-1:1 table
command -v glpsol >/dev/null || command -v scip >/dev/null \
  || { echo "No solver found: install GLPK (glpsol) in the jcvi environment"; exit 1; }

# remove any earlier 1:1 output so it is recomputed (the LAST alignment is kept and reused)
rm -rf "$RUN/work" "$RUN"/B73."$UNIT".lifted.1x1*

cd "$RUN"
$JCVI_PYTHON -m jcvi.compara.catalog ortholog B73 "$UNIT" \
  --cscore=.99 \
  --quota=1:1 \
  --no_strip_names \
  --no_dotplot \
  --cpus="${SLURM_CPUS_PER_TASK:-1}"

ls -l B73."$UNIT".*anchors*
python3 "$HERE/b73_table.py" --project "$PROJECT_DIR" "$UNIT"
echo "done $(date)"
