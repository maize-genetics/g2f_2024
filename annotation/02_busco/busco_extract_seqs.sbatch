#!/bin/bash
#SBATCH --job-name=busco_extract
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Step 02a: protein and CDS FASTA files from the whole-genome Helixer GFF3 of each unit (gffread),
# with checks. Sequence IDs are Helixer's own transcript IDs; nothing is renamed.
#
# Reads:  PROJECT_DIR/units.tsv, PROJECT_DIR/inputs/<unit>.scaffolds.fa,
#         PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.raw.gff3        (from step 01)
# Writes: PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.proteins.fa     terminal stop removed
#         PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.cds.fa          includes the stop codon
#         PROJECT_DIR/busco/extract_summary.tsv                           checks, one row per unit
#
# Usage (from the repository folder):
#   source config.sh && sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" 02_busco/busco_extract_seqs.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/config.sh"
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
GFFREAD=${GFFREAD:-gffread}

BASE=$PROJECT_DIR
mkdir -p "$BASE/busco"
SUMMARY=$BASE/busco/extract_summary.tsv

echo "gffread $($GFFREAD --version 2>&1)"
printf "unit\tgenes\tproteins\tcds\tprot_starts_M_pct\tprot_with_internal_stop\tcds_not_multiple_of_3\tmedian_protein_aa\tcheck\n" > "$SUMMARY"

while IFS=$'\t' read -r UNIT _; do
  [[ -n $UNIT ]] || continue
  echo "== $UNIT  $(date)"
  GENOME=$BASE/inputs/${UNIT}.scaffolds.fa
  D=$BASE/raw/$UNIT
  GFF=$D/${UNIT}.scaffolds.helixer.raw.gff3
  PROT=$D/${UNIT}.scaffolds.helixer.proteins.fa
  CDS=$D/${UNIT}.scaffolds.helixer.cds.fa
  for f in "$GENOME" "$GFF"; do [[ -s $f ]] || { echo "missing $f"; exit 1; }; done

  if [[ -s $PROT && -s $CDS ]]; then
    echo "   already extracted; re-checking only"
  else
    # gffread writes the stop codon as '.'; it is stripped from protein ends below
    $GFFREAD "$GFF" -g "$GENOME" -y "$PROT.tmp" -x "$CDS.tmp"
    python3 - "$PROT.tmp" "$PROT" <<'PY'
import sys
src, dst = sys.argv[1:3]
def recs(p):
    name, seq = None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield name, "".join(seq)
with open(dst, "w") as o:
    for name, s in recs(src):
        s = s.rstrip(".*")
        o.write(name + "\n")
        for i in range(0, len(s), 60):
            o.write(s[i:i+60] + "\n")
PY
    rm -f "$PROT.tmp"; mv "$CDS.tmp" "$CDS"
  fi

  # checks: one protein and one CDS per gene, same IDs, start/stop and frame statistics
  python3 - "$UNIT" "$GFF" "$PROT" "$CDS" >> "$SUMMARY" <<'PY'
import sys, statistics
unit, gff, prot, cds = sys.argv[1:5]
genes = 0
for l in open(gff):
    if l.startswith("#"): continue
    p = l.split("\t")
    if len(p) > 2 and p[2] == "gene": genes += 1
def recs(p):
    name, seq = None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name: yield name, "".join(seq)
P = dict(recs(prot)); C = dict(recs(cds))
startM = sum(1 for s in P.values() if s.startswith("M"))
internal = sum(1 for s in P.values() if "." in s or "*" in s)
not3 = sum(1 for s in C.values() if len(s) % 3)
med = statistics.median(len(s) for s in P.values()) if P else 0
ok = (len(P) == genes == len(C)) and set(P) == set(C)
print("\t".join(map(str, [unit, genes, len(P), len(C), f"{100*startM/max(len(P),1):.1f}",
                          internal, not3, f"{med:.0f}", "OK" if ok else "FAIL: counts or IDs differ"])))
PY
done < "$BASE/units.tsv"

echo
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
if grep -q "FAIL" "$SUMMARY"; then echo "Some units FAILED; see $SUMMARY"; exit 1; fi
echo "Done $(date)"
