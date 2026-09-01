#!/bin/bash
#SBATCH --job-name=merged_collage_pdfs
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/merged_collage_pdfs_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/merged_collage_pdfs_%j.err

# Build the four merged multi-page two-clone collages -- one PDF per SCNA
# holding every (sample x resolution) page in order -- and sync them into the
# existing two_clone_scna_collage_filtered folder on Drive.
#
# New targets (R/pipeline_targets_integration.R, inside the same tar_map as the
# per-resolution branches): two_clone_res_collages_filtered_merged_{1q,2p,6p,16q}
#   -> results/two_clone_scna_collage_filtered_merged_<scna>.pdf
# Pure qpdf concatenation of PDFs that already exist from job 11410385, so this
# is minutes of work; the bulk of the wall clock is hypoxia_rebuilt_gdrive's
# rclone pass. No numbatHelpers change, so no install_local step.
#
# One tar_make() at a time against _targets_r431 (AGENTS.md) -- check the queue
# before submitting:  squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"

module load r/4.4.1
module load rclone/1.74.3   # needed by hypoxia_rebuilt_gdrive
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Fail before the work rather than after it.
command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote configured"; exit 1; }

# The merged targets carry cue="always", so they need no invalidation; the
# per-resolution branches they read are current from job 11410385 and are
# deliberately NOT invalidated here -- re-rendering them is a 3 h job and
# nothing upstream of them changed.
Rscript run_targets_bg.R \
  two_clone_res_collages_filtered_merged_1q \
  two_clone_res_collages_filtered_merged_2p \
  two_clone_res_collages_filtered_merged_6p \
  two_clone_res_collages_filtered_merged_16q \
  hypoxia_rebuilt_gdrive

for f in results/two_clone_scna_collage_filtered_merged_*.pdf; do
  [ -f "$f" ] && echo "MERGED: $f ($(stat -c%s "$f") bytes)"
done

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"merged_collage_pdfs ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
