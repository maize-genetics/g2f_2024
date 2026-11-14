#!/bin/bash

#SBATCH --job-name=te_edta_helitron
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=5-00:00:00

# Step 06b (Helitron): detect Helitron transposable elements with EDTA_raw.pl.
# Runs at the same time as the other four raw jobs; all five outputs feed into te_edta_combine.
#
# Uses EDTA 2.3.0.
#
# Inputs:
#   $GENOME                         genome assembly FASTA
#
# Outputs:
#   edta_raw/$NAME.helitron.tar Helitron raw TE files  -> te_edta_combine

set -euo pipefail

# ── USER SETTINGS ────────────────────────────────────────────────────────────
# Edit these, or pass them at submission time, e.g.
#   sbatch --export=ALL,GENOME=/path/genome.fa,EDTA230_SIF=/path/edta.230.sif te_edta_helitron.sbatch
# Use the same GENOME and NAME for every TE job (te_submit.sh does this).

GENOME="${GENOME:-/path/to/genome.fa}"              # genome assembly FASTA
NAME="${NAME:-$(basename "${GENOME%.*}")}"          # sample name used for output files
SPECIES="${SPECIES:-Maize}"                         # EDTA --species: Maize, Rice, or others
EDTA230_SIF="${EDTA230_SIF:-}"                     # EDTA 2.3.0 image (docker://quay.io/biocontainers/edta:2.3.0--hdfd78af_0)
WORKDIR="${WORKDIR:-$PWD}"                          # where outputs are written (default: submission dir)
BIND="${BIND:-}"                                    # extra apptainer bind paths, comma-separated (optional)

# ── setup ────────────────────────────────────────────────────────────────────

type="helitron"
rawDir="$WORKDIR/edta_raw"
thr="${SLURM_CPUS_PER_TASK:-8}"

# Scratch space: a private folder on the node-local disk ($TMPDIR, or /tmp).
# /tmp is shared by all jobs on a node, so each job gets its own folder; it is
# deleted when the job ends (set KEEP_SCRATCH=1 to keep it for debugging).
scratch=$(mktemp -d "${TMPDIR:-/tmp}/edta_${NAME}_${type}.XXXXXX")
cleanup() { cd /; [[ "${KEEP_SCRATCH:-0}" == 1 ]] || rm -rf "$scratch"; }
trap cleanup EXIT

# Run a command inside an apptainer image, binding the work and scratch dirs
appt() {
    apptainer exec --cleanenv --bind "$WORKDIR,$scratch${BIND:+,$BIND}" "$@"
}
EDTA_raw() { appt "$EDTA230_SIF" EDTA_raw.pl "$@"; }

echo "Date: $(date --iso-8601=seconds)"
echo "Job ID: ${SLURM_JOB_ID:-NA}   Host: $HOSTNAME"
echo "GENOME=$GENOME"
echo "NAME=$NAME"
echo "SPECIES=$SPECIES"
echo "EDTA230_SIF=$EDTA230_SIF"
echo "WORKDIR=$WORKDIR"
echo "scratch=$scratch"
echo "threads=$thr"
echo

# ── input validation ─────────────────────────────────────────────────────────

command -v apptainer >/dev/null || { echo "Failed. apptainer not found (try: module load apptainer)"; exit 1; }
[[ -f "$GENOME" ]]      || { echo "Failed. Missing genome: $GENOME"; exit 1; }
[[ -f "$EDTA230_SIF" ]] || { echo "Failed. Missing EDTA 2.3.0 image: $EDTA230_SIF"; exit 1; }

mkdir -p "$rawDir"

# ── copy genome to scratch ───────────────────────────────────────────────────

genome_base=$(basename "$GENOME")
cp "$GENOME" "$scratch/$genome_base"

# ── EDTA_raw ─────────────────────────────────────────────────────────────────

echo "Start EDTA_raw --type $type"

cd "$scratch"
unset -f which    # some module systems define a 'which' function that breaks EDTA

EDTA_raw \
    --genome  "$scratch/$genome_base" \
    --species "$SPECIES" \
    --type    "$type" \
    --threads "$thr"

echo "EDTA_raw $type complete"
echo

# ── validate and tar output ──────────────────────────────────────────────────
# EDTA names its outputs <genome>.mod*; te_edta_combine needs all of them.

tarout="$rawDir/$NAME.$type.tar"

shopt -s nullglob
modfiles=("$scratch/${genome_base}.mod"*)
[[ ${#modfiles[@]} -gt 0 ]] || { echo "Failed. No ${genome_base}.mod* files found after EDTA_raw $type"; exit 1; }

tar -C "$scratch" -cf "$tarout" "${modfiles[@]##*/}"

echo "Helitron tar: $tarout"
echo "$(date --iso-8601=seconds) done"
