#!/bin/bash
#SBATCH --job-name=force_rclone_hyp
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/force_rclone_hyp_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/force_rclone_hyp_%j.err

# hypoxia_rebuilt_gdrive was skipped after the collage/summary rebuild because its
# dependencies (hypoxia_summaries, annotated_heatmap_collages) return the SAME
# deterministic results/*.pdf path strings each run -- targets hashes the string,
# not the file content, so the new-label PDFs on disk didn't change the dep hashes.
# Force a re-sync by invalidating hypoxia_rebuilt_gdrive, then tar_make it: it
# re-reads the updated files from disk and rclones them to Google Drive.

module load r/4.4.1 rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(any_of("hypoxia_rebuilt_gdrive"))
cat("invalidated hypoxia_rebuilt_gdrive\n")
' || { echo "FATAL: invalidation errored"; exit 1; }

Rscript run_targets_bg.R hypoxia_rebuilt_gdrive
