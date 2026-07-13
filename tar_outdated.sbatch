#!/bin/bash
#SBATCH --job-name=tar_outdated
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/tar_outdated_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/tar_outdated_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
  suppressPackageStartupMessages({library(targets); library(numbatHelpers)})
  od <- targets::tar_outdated(callr_function = NULL)
  cat("=== OUTDATED TARGETS (", length(od), ") ===\n", sep = "")
  cat(sort(od), sep = "\n")
  cat("\n=== is cluster_dictionary outdated? ===\n")
  cat("cluster_dictionary" %in% od, "\n")
  cat("=== is sample_summaries outdated? ===\n")
  cat("sample_summaries" %in% od, "\n")
'
