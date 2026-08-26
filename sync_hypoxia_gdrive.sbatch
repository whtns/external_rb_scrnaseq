#!/bin/bash
#SBATCH --job-name=sync_hypoxia_gdrive
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/sync_hypoxia_gdrive_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/sync_hypoxia_gdrive_%j.err

# Re-run hypoxia_rebuilt_gdrive only. The target is already up to date from the
# last build, so it must be invalidated to force the sync; nothing upstream is
# touched, so no figure is rebuilt.
#
# One tar_make() at a time against _targets_r431 (AGENTS.md).

module load r/4.4.1
module load rclone/1.74.3
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
tar_invalidate(any_of("hypoxia_rebuilt_gdrive"))
cat("invalidated hypoxia_rebuilt_gdrive\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

Rscript run_targets_bg.R hypoxia_rebuilt_gdrive

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"sync_hypoxia_gdrive ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
