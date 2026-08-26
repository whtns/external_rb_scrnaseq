#!/bin/bash
#SBATCH --job-name=maj_rds
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --array=0-26%8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/maj_rds_%A_%a.log

# Rebuild <sample>_numbat.rds at the round selected by the MAJORITY-OF-ROUNDS
# rule (see docs/numbat_round_policy.md sections 4b/4c).
#
# An RB event is reported when a majority of a sample's consensus rounds call
# it; the RDS is then built from the round whose own event set equals that
# majority set, so the object displays exactly what the table reports. Such a
# round exists for all 36 scored samples.
#
# Sample list is read from results/rds_rebuild_list_20260826.txt, which
# src/write_majority_manifest.R writes: the samples whose NEW round differs
# from the round their on-disk RDS was built at. Reading the list from disk
# (rather than deriving it) keeps the array index -> sample mapping fixed even
# if the manifest is regenerated mid-flight.
#
# SCOPE, per AGENTS.md: SRX only, never SRR, and never the three held-out
# samples. The list contains none of them; the guard below is belt-and-braces.
#
# Existing RDS were backed up to output/numbat_rds_backup_20260826/ before this
# job was first run. The previous manifest is at
# results/numbat_selected_round_premajority_20260826.csv.
#
# Uses the production container (R 4.3.1, numbat 1.5.2) to match the
# _targets_r431 store, not the module R stack.

set -euo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
SRC_DIR=$PROJ/output/numbat_sridhar
MANIFEST=$PROJ/results/numbat_selected_round.csv
# Override with --export=ALL,RDS_LIST=<path> to rebuild a different subset
# (e.g. results/rds_rebuild_list_pass2.txt after a verification pass). Set the
# --array range to match the list length.
LIST=${RDS_LIST:-$PROJ/results/rds_rebuild_list_20260826.txt}

# Held out at the user's request; see AGENTS.md.
EXCLUDE="SRX11133592 SRX11133593 SRX11133594"

cd "$PROJ"

mapfile -t SAMPLES < "$LIST"

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#SAMPLES[@]}" ]; then
  echo "array index ${SLURM_ARRAY_TASK_ID} beyond sample list (${#SAMPLES[@]}); nothing to do"
  exit 0
fi

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

case " $EXCLUDE " in
  *" $S "*) echo "REFUSING: $S is held out per AGENTS.md"; exit 1 ;;
esac
case "$S" in
  SRX*) ;;
  *) echo "REFUSING: $S is not an SRX sample"; exit 1 ;;
esac

K=$(awk -F, -v s="$S" 'NR==1{for(i=1;i<=NF;i++){if($i=="sample_id")c=i;if($i=="selected_round")r=i};next} $c==s{print $r; exit}' "$MANIFEST")
PREV=$(awk -F, -v s="$S" 'NR==1{for(i=1;i<=NF;i++){if($i=="sample_id")c=i;if($i=="prev_selected_round")r=i};next} $c==s{print $r; exit}' "$MANIFEST")

echo "=== $S : rebuilding _numbat.rds  round ${PREV:-?} -> ${K:-<not in manifest>} ==="
echo "target: $SRC_DIR/${S}_numbat.rds"
echo "backup: $PROJ/output/numbat_rds_backup_20260826/${S}_numbat.rds"

if [ ! -f "$PROJ/output/numbat_rds_backup_20260826/${S}_numbat.rds" ]; then
  echo "REFUSING: no backup on disk for $S"
  exit 1
fi

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
