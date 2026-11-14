#!/bin/bash

#SBATCH --job-name=te_repeatmasker
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=5-00:00:00

# Step 06d: run RepeatMasker against the EDTA TE library (from te_edta_combine) and make
# a soft-masked genome FASTA. Run after te_edta_combine has finished.
#
# Masking: -q -div 40 -cutoff 225 balances speed and sensitivity.
# The RepeatMasker GFF (not EDTA's annotation) drives bedtools maskfasta.
# bedtools is taken from the EDTA 2.3.0 image by default (set BEDTOOLS_SIF to
# use a different image).
#
# Inputs:
#   $GENOME                             genome assembly FASTA (same as the raw jobs)
#   <genome>.mod.EDTA.TElib.fa          EDTA TE library (te_edta_combine)
#
# Outputs:
#   repeatmasker/$NAME.tar           RepeatMasker output files (.out, .gff, .tbl, ...)
#   $NAME.softmasked.fa                 soft-masked genome FASTA -> $WORKDIR

set -euo pipefail

# ── USER SETTINGS ────────────────────────────────────────────────────────────
# Edit these, or pass them at submission time, e.g.
#   sbatch --export=ALL,GENOME=/path/genome.fa,RM_SIF=/path/repeatmasker.sif,EDTA230_SIF=/path/edta.230.sif te_repeatmasker.sbatch
# Use the same GENOME and NAME for every TE job (te_submit.sh does this).

GENOME="${GENOME:-/path/to/genome.fa}"              # genome assembly FASTA
NAME="${NAME:-$(basename "${GENOME%.*}")}"          # sample name used for output files
RM_SIF="${RM_SIF:-}"                               # RepeatMasker 4.2.3 image (docker://quay.io/biocontainers/repeatmasker:4.2.3--pl5321hdfd78af_0)
EDTA230_SIF="${EDTA230_SIF:-}"                     # EDTA 2.3.0 image, used only for bedtools
BEDTOOLS_SIF="${BEDTOOLS_SIF:-$EDTA230_SIF}"        # any image with bedtools; defaults to the EDTA image
WORKDIR="${WORKDIR:-$PWD}"                          # where outputs are written (default: submission dir)
BIND="${BIND:-}"                                    # extra apptainer bind paths, comma-separated (optional)

# ── setup ────────────────────────────────────────────────────────────────────

rmaskDir="$WORKDIR/repeatmasker"
genome_base=$(basename "$GENOME")
lib="$WORKDIR/${genome_base}.mod.EDTA.TElib.fa"     # written by te_edta_combine
softmasked="$WORKDIR/$NAME.softmasked.fa"
thr="${SLURM_CPUS_PER_TASK:-8}"

# Scratch space: a private folder on the node-local disk ($TMPDIR, or /tmp).
# /tmp is shared by all jobs on a node, so each job gets its own folder; it is
# deleted when the job ends (set KEEP_SCRATCH=1 to keep it for debugging).
scratch=$(mktemp -d "${TMPDIR:-/tmp}/edta_${NAME}_rmask.XXXXXX")
cleanup() { cd /; [[ "${KEEP_SCRATCH:-0}" == 1 ]] || rm -rf "$scratch"; }
trap cleanup EXIT

# Run a command inside an apptainer image, binding the work and scratch dirs
appt() {
    apptainer exec --cleanenv --bind "$WORKDIR,$scratch${BIND:+,$BIND}" "$@"
}
RepeatMasker() { appt "$RM_SIF" RepeatMasker "$@"; }
bedtools()     { appt "$BEDTOOLS_SIF" bedtools "$@"; }

echo "Date: $(date --iso-8601=seconds)"
echo "Job ID: ${SLURM_JOB_ID:-NA}   Host: $HOSTNAME"
echo "GENOME=$GENOME"
echo "NAME=$NAME"
echo "lib=$lib"
echo "RM_SIF=$RM_SIF"
echo "BEDTOOLS_SIF=$BEDTOOLS_SIF"
echo "WORKDIR=$WORKDIR"
echo "scratch=$scratch"
echo "threads=$thr"
echo

# ── input validation ─────────────────────────────────────────────────────────

command -v apptainer >/dev/null || { echo "Failed. apptainer not found (try: module load apptainer)"; exit 1; }
[[ -f "$GENOME" ]]      || { echo "Failed. Missing genome: $GENOME"; exit 1; }
[[ -f "$lib" ]]         || { echo "Failed. Missing EDTA TE library: $lib (run te_edta_combine first)"; exit 1; }
[[ -f "$RM_SIF" ]]      || { echo "Failed. Missing RepeatMasker image: $RM_SIF"; exit 1; }
[[ -f "$BEDTOOLS_SIF" ]] || { echo "Failed. Missing bedtools image: $BEDTOOLS_SIF"; exit 1; }

mkdir -p "$rmaskDir"

# ── copy genome to scratch ───────────────────────────────────────────────────

cp "$GENOME" "$scratch/$genome_base"
cd "$scratch"

# ── RepeatMasker ─────────────────────────────────────────────────────────────

echo "Start RepeatMasker"

unset -f which    # some module systems define a 'which' function that breaks RepeatMasker
RepeatMasker \
    -e      ncbi \
    -pa     "$thr" \
    -q \
    -div    40 \
    -cutoff 225 \
    -gff \
    -lib    "$lib" \
    "$scratch/$genome_base"

# Fix missing TE subclass: RepeatMasker emits a bare "DNA" classification;
# downstream tools (e.g. LAI) expect "DNA/unknown".
perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' "$scratch/${genome_base}.out"

echo "RepeatMasker complete"
echo

# ── archive RepeatMasker outputs ─────────────────────────────────────────────

tar -C "$scratch" -cf "$rmaskDir/$NAME.tar" \
    "${genome_base}.out" \
    "${genome_base}.out.gff" \
    "${genome_base}.tbl" \
    "${genome_base}.masked" \
    "${genome_base}.cat.gz"

echo "RepeatMasker tar: $rmaskDir/$NAME.tar"
echo

# ── soft-masked FASTA ────────────────────────────────────────────────────────
# Lower-case every base covered by a RepeatMasker hit.

echo "Generating soft-masked FASTA with bedtools maskfasta..."

bedtools maskfasta \
    -soft \
    -fi  "$scratch/$genome_base" \
    -bed "$scratch/${genome_base}.out.gff" \
    -fo  "$softmasked"

echo "Soft-masked FASTA: $softmasked"
echo "$(date --iso-8601=seconds) done"
