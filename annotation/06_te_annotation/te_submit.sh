#!/bin/bash
# Step 06: submit the TE annotation of one or more genomes as a chain of Slurm jobs (each waits
# for the one before):
#
#   te_edta_ltr      \
#   te_edta_tir       |
#   te_edta_helitron  +--> te_edta_combine --> te_repeatmasker --> te_restore_ids
#   te_edta_line      |
#   te_edta_sine     /
#
# Each genome runs in PROJECT_DIR/te_annotation/runs/<unit>/ (job logs in its logs/ folder).
#
# Usage (from the repository folder, after te_prep_inputs.py):
#   bash 06_te_annotation/te_submit.sh                    # every unit in units.tsv
#   bash 06_te_annotation/te_submit.sh UNIT1 UNIT2        # some units
#   START=combine bash 06_te_annotation/te_submit.sh UNIT1        # resume: combine, repeatmasker or restore
set -euo pipefail

source ./config.sh
: "${PROJECT_DIR:?set PROJECT_DIR in config.sh}"
SCRIPTS=$(cd "$(dirname "$0")" && pwd)
TE=$PROJECT_DIR/te_annotation

SPECIES=${SPECIES:-Maize}        # EDTA --species
SENSITIVE=${SENSITIVE:-1}        # 1 = EDTA --sensitive 1 (adds RepeatModeler; slower)
START=${START:-raw}
case $START in raw|combine|repeatmasker|restore) ;; *) echo "START must be raw, combine, repeatmasker or restore"; exit 1 ;; esac

RM_SIF=$REPEATMASKER_SIF
for v in EDTA230_SIF EDTA222_SIF RM_SIF; do
  [[ -f ${!v} ]] || { echo "Set $v in config.sh (now: '${!v}')"; exit 1; }
done

# account, partition and e-mail from config.sh, if set
SB=(sbatch --parsable)
[[ -n ${SLURM_ACCOUNT:-} ]] && SB+=(-A "$SLURM_ACCOUNT")
[[ -n ${SLURM_PARTITION:-} ]] && SB+=(-p "$SLURM_PARTITION")
[[ -n ${EMAIL:-} ]] && SB+=(--mail-type=FAIL,END --mail-user="$EMAIL")

UNITS=("$@")
[[ ${#UNITS[@]} -gt 0 ]] || mapfile -t UNITS < <(cut -f1 "$PROJECT_DIR/units.tsv" | sed '/^$/d')

for UNIT in "${UNITS[@]}"; do
  GENOME=$TE/inputs/$UNIT.nuclear.fa
  CDS=$PROJECT_DIR/nuclear_genome/$UNIT/$UNIT.nuclear.helixer.cds.fa   # used by EDTA to remove gene sequences
  WORKDIR=$TE/runs/$UNIT
  NAME=$UNIT
  [[ -s $GENOME ]] || { echo "$UNIT: missing $GENOME (run te_prep_inputs.py)"; exit 1; }
  [[ -s $CDS ]]    || { echo "$UNIT: missing $CDS (run step 04)"; exit 1; }
  mkdir -p "$WORKDIR/logs"
  export GENOME NAME SPECIES SENSITIVE CDS EDTA222_SIF EDTA230_SIF RM_SIF WORKDIR TE SCRIPTS

  submit() {   # submit <script> [job ids to wait for ...]
    local script=$1; shift
    local dep=()
    [[ $# -gt 0 ]] && dep=(--dependency=afterok:$(IFS=:; echo "$*"))
    "${SB[@]}" ${dep[@]+"${dep[@]}"} --chdir "$WORKDIR" --job-name "${NAME}_${script%.sbatch}" \
      --export=ALL "$SCRIPTS/$script"
  }

  echo "== $UNIT   genome $GENOME   out $WORKDIR"
  raw=(); comb=(); rmj=()
  if [[ $START == raw ]]; then
    for s in te_edta_ltr te_edta_tir te_edta_helitron te_edta_line te_edta_sine; do
      j=$(submit "$s.sbatch"); echo "   $s $j"; raw+=("$j")
    done
  fi
  if [[ $START == raw || $START == combine ]]; then
    j=$(submit te_edta_combine.sbatch ${raw[@]+"${raw[@]}"}); echo "   te_edta_combine $j"; comb=("$j")
  fi
  if [[ $START != restore ]]; then
    j=$(submit te_repeatmasker.sbatch ${comb[@]+"${comb[@]}"}); echo "   te_repeatmasker $j"; rmj=("$j")
  fi
  j=$(submit te_restore_ids.sbatch ${rmj[@]+"${rmj[@]}"}); echo "   te_restore_ids $j"
done
