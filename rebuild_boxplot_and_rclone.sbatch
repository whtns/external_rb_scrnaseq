#!/bin/bash
#SBATCH --job-name=boxplot_rclone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/boxplot_rclone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/boxplot_rclone_%j.err

# Regenerate the per-cluster mean_score boxplots from the current
# hypoxia_split_log_all.csv, then rclone them to Google Drive.
#
# The boxplot script is a standalone Rscript (not a target), so run it first;
# hypoxia_rebuilt_gdrive then picks the PDF up via its mean_score_boxplot entry.
# Only hypoxia_rebuilt_gdrive is invalidated -- the split itself does not re-run.

module load r/4.4.1 rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/plot_hypoxia_mean_score_boxplots.R || { echo "FATAL: boxplot script failed"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(any_of("hypoxia_rebuilt_gdrive"))
cat("invalidated hypoxia_rebuilt_gdrive\n")
' || { echo "FATAL: invalidation errored"; exit 1; }

Rscript run_targets_bg.R hypoxia_rebuilt_gdrive
