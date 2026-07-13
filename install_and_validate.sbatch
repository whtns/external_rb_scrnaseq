#!/bin/bash
#SBATCH --job-name=install_validate
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/install_validate_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/install_validate_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force = TRUE, upgrade = "never")' && \
Rscript -e '
  suppressPackageStartupMessages(library(numbatHelpers))
  cat("build_cluster_dictionary exported:", exists("build_cluster_dictionary"), "\n")
  cat("classify_cluster_abbreviation exported:", exists("classify_cluster_abbreviation"), "\n")
  library(targets)
  targets::tar_validate(callr_function = NULL)
  cat("tar_validate: OK (no errors / no cycle)\n")
'
