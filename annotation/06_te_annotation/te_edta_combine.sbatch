#!/bin/bash

#SBATCH --job-name=te_edta_combine
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=5-00:00:00

# Step 06c: combine the per-type EDTA_raw results and finish the EDTA pipeline.
# Run after all five raw detection jobs have finished.
#
# Extracts the five per-type tars into scratch, then runs EDTA.pl with
# --anno 0 --evaluate 0 so it skips raw detection and runs filter + finish
# only, using the pre-computed outputs.
#
# Uses EDTA 2.3.0. RepeatMasker and soft-masking are done in te_repeatmasker.
#
# Inputs:
#   $GENOME                              genome assembly FASTA (same as the raw jobs)
#   edta_raw/$NAME.{type}.tar         per-type raw TE files (the five raw jobs)
#   $CDS                                 CDS FASTA to filter gene TEs (optional)
#
# Outputs:
#   <genome>.mod.EDTA.TElib.fa           TE consensus library   -> $WORKDIR (te_repeatmasker)
#   <genome>.mod.EDTA.intact.gff3        intact TE annotation   -> $WORKDIR
#   edta/$NAME.tar                    all EDTA .mod* files   (archived)

set -euo pipefail

# ── USER SETTINGS ────────────────────────────────────────────────────────────
# Edit these, or pass them at submission time, e.g.
#   sbatch --export=ALL,GENOME=/path/genome.fa,EDTA230_SIF=/path/edta.230.sif te_edta_combine.sbatch
# Use the same GENOME and NAME for every TE job (te_submit.sh does this).

GENOME="${GENOME:-/path/to/genome.fa}"              # genome assembly FASTA
NAME="${NAME:-$(basename "${GENOME%.*}")}"          # sample name used for output files
SPECIES="${SPECIES:-Maize}"                         # EDTA --species: Maize, Rice, or others
SENSITIVE="${SENSITIVE:-1}"                         # 1 = sensitive (thorough, slow); 0 = fast
CDS="${CDS:-}"                                      # optional CDS FASTA to reduce false TE calls in genes
EDTA230_SIF="${EDTA230_SIF:-}"                     # EDTA 2.3.0 image (docker://quay.io/biocontainers/edta:2.3.0--hdfd78af_0)
WORKDIR="${WORKDIR:-$PWD}"                          # where outputs are written (default: submission dir)
BIND="${BIND:-}"                                    # extra apptainer bind paths, comma-separated (optional)

# ── setup ────────────────────────────────────────────────────────────────────

rawDir="$WORKDIR/edta_raw"
edtaDir="$WORKDIR/edta"
thr="${SLURM_CPUS_PER_TASK:-8}"

# Scratch space: a private folder on the node-local disk ($TMPDIR, or /tmp).
# /tmp is shared by all jobs on a node, so each job gets its own folder; it is
# deleted when the job ends (set KEEP_SCRATCH=1 to keep it for debugging).
scratch=$(mktemp -d "${TMPDIR:-/tmp}/edta_${NAME}_combine.XXXXXX")
cleanup() { cd /; [[ "${KEEP_SCRATCH:-0}" == 1 ]] || rm -rf "$scratch"; }
trap cleanup EXIT

# Run a command inside an apptainer image, binding the work and scratch dirs
appt() {
    apptainer exec --cleanenv --bind "$WORKDIR,$scratch${BIND:+,$BIND}" "$@"
}
EDTA.pl() { appt "$EDTA230_SIF" EDTA.pl "$@"; }

echo "Date: $(date --iso-8601=seconds)"
echo "Job ID: ${SLURM_JOB_ID:-NA}   Host: $HOSTNAME"
echo "GENOME=$GENOME"
echo "NAME=$NAME"
echo "SPECIES=$SPECIES"
echo "SENSITIVE=$SENSITIVE"
echo "CDS=${CDS:-<not set>}"
echo "EDTA230_SIF=$EDTA230_SIF"
echo "WORKDIR=$WORKDIR"
echo "scratch=$scratch"
echo "threads=$thr"
echo

# ── input validation ─────────────────────────────────────────────────────────

command -v apptainer >/dev/null || { echo "Failed. apptainer not found (try: module load apptainer)"; exit 1; }
[[ -f "$GENOME" ]]      || { echo "Failed. Missing genome: $GENOME"; exit 1; }
[[ -f "$EDTA230_SIF" ]] || { echo "Failed. Missing EDTA 2.3.0 image: $EDTA230_SIF"; exit 1; }

for type in ltr tir helitron line sine; do
    tar="$rawDir/$NAME.$type.tar"
    [[ -f "$tar" ]] || { echo "Failed. Missing raw tar: $tar (run the five raw jobs first)"; exit 1; }
done

if [[ -n "$CDS" ]]; then
    [[ -f "$CDS" ]] || { echo "Failed. CDS is set but not found: $CDS"; exit 1; }
fi

mkdir -p "$edtaDir"

# ── copy genome and extract raw tars to scratch ──────────────────────────────

genome_base=$(basename "$GENOME")
cp "$GENOME" "$scratch/$genome_base"
[[ -n "$CDS" ]] && cp "$CDS" "$scratch/cds.fa"

echo "Extracting per-type raw tars to scratch..."
for type in ltr tir helitron line sine; do
    echo "  $NAME.$type.tar"
    tar -C "$scratch" -xf "$rawDir/$NAME.$type.tar"
done
echo

# ── EDTA combine ─────────────────────────────────────────────────────────────

echo "Start EDTA.pl (filter + finish, using pre-computed raw outputs)"

cd "$scratch"
unset -f which    # some module systems define a 'which' function that breaks EDTA

edta_cmd=(
    EDTA.pl
    --genome    "$scratch/$genome_base"
    --species   "$SPECIES"
    --sensitive "$SENSITIVE"
    --anno      0
    --evaluate  0
    --overwrite 0
    --threads   "$thr"
)
[[ -n "$CDS" ]] && edta_cmd+=(--cds "$scratch/cds.fa")

echo "Running: ${edta_cmd[*]}"
"${edta_cmd[@]}"

echo "EDTA.pl combine complete"
echo

# ── validate, copy, and tar EDTA outputs ─────────────────────────────────────

lib="${genome_base}.mod.EDTA.TElib.fa"
anno="${genome_base}.mod.EDTA.intact.gff3"

[[ -f "$scratch/$lib" ]]  || { echo "Failed. EDTA TE library not found: $scratch/$lib"; exit 1; }
[[ -f "$scratch/$anno" ]] || { echo "Failed. EDTA intact annotation not found: $scratch/$anno"; exit 1; }

cp "$scratch/$lib"  "$WORKDIR/$lib"
cp "$scratch/$anno" "$WORKDIR/$anno"

echo "TE library:         $WORKDIR/$lib"
echo "TE annotation GFF3: $WORKDIR/$anno"
echo "Annotated features: $(grep -vc '^#' "$WORKDIR/$anno" || true)"
echo

shopt -s nullglob
modfiles=("$scratch/${genome_base}.mod"*)
tar -C "$scratch" -cf "$edtaDir/$NAME.tar" "${modfiles[@]##*/}"

echo "EDTA combine tar: $edtaDir/$NAME.tar"
echo "$(date --iso-8601=seconds) done"
