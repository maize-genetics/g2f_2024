#!/bin/bash
#SBATCH --job-name=helixer_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 01a: prepare the genome assemblies for Helixer (whole genomes, organelle sequence kept).
#
# Reads ASSEMBLIES_TSV (see assemblies_example.tsv): one row per genome assembly ("unit") with
#   unit, source label, FASTA path (.fa or .fa.gz).
# Each FASTA is copied (or decompressed) to PROJECT_DIR/inputs/<unit>.scaffolds.fa; its number of
# sequences, total bp and sha256 are recorded, and an uncompressed copy is checked against its source.
# The source files are never written.
#
# Writes: PROJECT_DIR/units.tsv                    unit -> input FASTA (read by every later step)
#         PROJECT_DIR/inputs/<unit>.scaffolds.fa
#         PROJECT_DIR/inputs/prep_check.tsv         per unit: sequences, bp, sha256, status OK / FAIL
#
# Usage (from the repository folder):
#   source config.sh && mkdir -p logs
#   sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" 01_helixer/helixer_prep_inputs.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}" "${ASSEMBLIES_TSV:?set ASSEMBLIES_TSV in config.sh}"

BASE=$PROJECT_DIR
mkdir -p "$BASE"/{inputs,raw,final,qc,logs}
CHECK=$BASE/inputs/prep_check.tsv
printf "unit\tsource\tstatus\trecords\tbp\tsha256\n" > "$CHECK"
: > "$BASE/units.tsv"

count() {  # records and bp of a FASTA
  awk '/^>/{n++; next} {bp+=length($0)} END{print n+0"\t"bp+0}' "$1"
}

while IFS=$'\t' read -r UNIT SOURCE SRC _; do
  [[ -z $UNIT || $UNIT == unit || $UNIT == \#* ]] && continue
  echo "== $UNIT"
  [[ -s $SRC ]] || { echo "   missing FASTA for $UNIT: $SRC"; exit 1; }

  OUTFA=$BASE/inputs/$UNIT.scaffolds.fa
  if [[ ! -s $OUTFA ]]; then
    if [[ $SRC == *.gz ]]; then zcat "$SRC" > "$OUTFA.tmp"; else cp "$SRC" "$OUTFA.tmp"; fi
    mv "$OUTFA.tmp" "$OUTFA"
  fi
  read -r GOT_REC GOT_BP < <(count "$OUTFA")

  STATUS=OK
  [[ $GOT_REC -gt 0 && $GOT_BP -gt 0 ]] || STATUS=FAIL
  [[ $SRC != *.gz ]] && ! cmp -s "$SRC" "$OUTFA" && STATUS=FAIL
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$UNIT" "$SOURCE" "$STATUS" "$GOT_REC" "$GOT_BP" \
    "$(sha256sum "$OUTFA" | cut -d' ' -f1)" >> "$CHECK"
  printf "%s\t%s\n" "$UNIT" "$OUTFA" >> "$BASE/units.tsv"
done < "$ASSEMBLIES_TSV"

echo
column -t "$CHECK" 2>/dev/null || cat "$CHECK"
echo "units: $(wc -l < "$BASE/units.tsv")"
if grep -q $'\tFAIL\t' "$CHECK"; then echo "One or more units FAILED; see $CHECK"; exit 1; fi
echo "All inputs ready."
