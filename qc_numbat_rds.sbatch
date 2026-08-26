#!/bin/bash
#SBATCH --job-name=qc_numbat_rds
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/qc_numbat_rds_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/qc_numbat_rds_%j.err

# Build the numbat RDS diagnostics only:
#   numbat_rds_qc         one branch per *_numbat.rds -> summary row + per-sample PDF
#   numbat_rds_qc_table   results/numbat_rds_qc.csv
#   numbat_rds_qc_report  results/numbat_rds_qc_report.pdf (all pages combined)
#
# SRX samples only (39 of the 71 RDS). Upstream is numbat_rds_srx only, so this
# does NOT wake the Seurat/diffex/figure half. It still runs tar_make against
# _targets_r431, so the
# one-job rule applies: check `squeue -u $USER` for a running run_targets_bg.R or
# crew-worker-* before submitting.

set -euo pipefail

module load gcc/13.3.0 r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# The QC target calls numbatHelpers::qc_numbat_rds() inside crew workers, which
# are fresh R sessions on other nodes and load the INSTALLED package. load_all()
# in the controller is not enough, so install before running.
Rscript -e "devtools::install_local('/project2/cobrinik_1090/rpkgs/numbat_helpers', force = TRUE, upgrade = 'never')"

# targets hashes each target's COMMAND, not the internals of the packages that
# command calls. numbatHelpers::qc_numbat_rds() changing (e.g. the arm-coverage
# floor in numbat_rb_events) therefore leaves these targets looking up to date
# and they get skipped. Invalidate explicitly so a rerun really re-inspects.
# Cheap to do unconditionally -- the whole pattern is ~12 minutes.
Rscript -e "suppressPackageStartupMessages(library(targets)); \
  targets::tar_config_set(store = '_targets_r431'); \
  targets::tar_invalidate(targets::any_of(c('numbat_rds_qc','numbat_rds_qc_table','numbat_rds_qc_report'))); \
  cat('invalidated QC targets\n')"

Rscript run_targets_bg.R numbat_rds_srx numbat_rds_qc numbat_rds_qc_table numbat_rds_qc_report
