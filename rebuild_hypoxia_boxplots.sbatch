#!/bin/bash
#SBATCH --job-name=hypoxia_boxplots
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/hypoxia_boxplots_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/hypoxia_boxplots_%j.err

# Refresh results/hypoxia_cluster_split/hypoxia_mean_score_boxplots.pdf from the
# CURRENT split logic. All 33 per-sample *_hypoxia_split_log.csv were rewritten
# 2026-07-17 by a full re-split, but the collated hypoxia_split_log_all.csv and the
# boxplot PDF are both from 2026-07-13 -- stale. Re-collate exactly as the
# hypoxia_split_log_collated target does (plain bind_rows of every per-sample log,
# no dedup/filter), then regenerate the plot. Standalone: touches only these two
# result files, not _targets_r431, so it is safe to run beside the collage rebuild.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# 1. Re-collate the per-sample logs (mirror of hypoxia_split_log_collated).
Rscript -e '
dir <- "results/hypoxia_cluster_split"
csvs <- list.files(dir, pattern = "_hypoxia_split_log\\.csv$", full.names = TRUE)
stopifnot(length(csvs) > 0)
out <- file.path(dir, "hypoxia_split_log_all.csv")
readr::write_csv(dplyr::bind_rows(lapply(csvs, readr::read_csv, show_col_types = FALSE)), out)
cat("re-collated", length(csvs), "per-sample logs ->", out, "\n")
' || { echo "FATAL: collation failed"; exit 1; }

# 2. Regenerate the boxplot PDF from the fresh CSV.
Rscript src/plot_hypoxia_mean_score_boxplots.R
