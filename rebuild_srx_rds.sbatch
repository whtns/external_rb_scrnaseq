#!/bin/bash
#SBATCH --job-name=srx_rds
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --array=0-35%8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/srx_rds_%A_%a.log

# Rebuild <sample>_numbat.rds for the SRX samples at their SELECTED consensus
# round (results/numbat_selected_round.csv) instead of numbat's i = 2 default.
#
# SCOPE, per AGENTS.md: SRX only, never SRR, and never the three held-out
# samples below. The list is derived here rather than passed in, so the array
# index -> sample mapping is reproducible from this script alone.
#
# Existing RDS files were backed up to output/numbat_rds_backup_20260820/
# before this job was first run.
#
# Uses the production container (R 4.3.1, numbat 1.5.2) to match the
# _targets_r431 store, not the module R stack.

set -euo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
SRC_DIR=$PROJ/output/numbat_sridhar
MANIFEST=$PROJ/results/numbat_selected_round.csv

# Held out at the user's request; see AGENTS.md.
EXCLUDE="SRX11133592 SRX11133593 SRX11133594"

cd "$PROJ"

# SRX sample dirs with a done.txt, minus the exclusions, alphabetical.
SAMPLES=()
for d in "$SRC_DIR"/SRX*/; do
  s=$(basename "$d")
  [ -f "$d/done.txt" ] || continue
  case " $EXCLUDE " in *" $s "*) continue ;; esac
  SAMPLES+=("$s")
done
IFS=$'\n' SAMPLES=($(sort <<<"${SAMPLES[*]}")); unset IFS

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#SAMPLES[@]}" ]; then
  echo "array index ${SLURM_ARRAY_TASK_ID} beyond sample list (${#SAMPLES[@]}); nothing to do"
  exit 0
fi

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
K=$(awk -F, -v s="$S" 'NR==1{for(i=1;i<=NF;i++){if($i=="sample_id")c=i;if($i=="selected_round")r=i};next} $c==s{print $r; exit}' "$MANIFEST")

echo "=== $S : rebuilding _numbat.rds at selected round ${K:-<not in manifest>} ==="
echo "target: $SRC_DIR/${S}_numbat.rds"

module purge
module load apptainer

APP=$(command -v apptainer)
SIF=$PROJ/pipeline/containers/numbat-pipeline.sif

export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

"$APP" exec --bind /project2 "$SIF" /usr/bin/Rscript \
  pipeline/scripts/process_numbat_rds.R \
  "done_file='$SRC_DIR/$S/done.txt'" \
  "selection_manifest='$MANIFEST'"

echo "=== $S done; rds now:"
ls -la "$SRC_DIR/${S}_numbat.rds"
