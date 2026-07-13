#!/bin/bash
#SBATCH --job-name=split_log_rclone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/split_log_rclone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/split_log_rclone_%j.err

# Re-collate results/hypoxia_cluster_split/hypoxia_split_log_all.csv from the
# per-sample *_hypoxia_split_log.csv files and push it to Google Drive.
#
# hypoxia_split_log_collated globs the per-sample logs and binds them, so
# invalidating it forces a fresh re-collation. hypoxia_partition_paths is NOT
# invalidated -- it stays cached (cue(command = FALSE)), so the expensive split
# itself does not re-run; only the collation + rclone do.

module load r/4.4.1 rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(any_of(c("hypoxia_split_log_collated", "hypoxia_rebuilt_gdrive")))
cat("invalidated hypoxia_split_log_collated + hypoxia_rebuilt_gdrive\n")
' || { echo "FATAL: invalidation errored"; exit 1; }

Rscript run_targets_bg.R hypoxia_rebuilt_gdrive
