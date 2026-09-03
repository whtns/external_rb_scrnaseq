#!/bin/bash
#SBATCH --job-name=smoke_bulk_panel
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_bulk_panel_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_bulk_panel_%j.err

# Smoke test for plot_numbat_bulk_clones() (github #41 / #42).
# Three samples, no targets store touched, writes only
# results/numbat_bulk_clones_SMOKE/.
#
# 64 GB: bulk_clones runs to millions of rows on the larger samples, which is
# why pipeline_targets_qc.R puts the equivalent work on the heavy controller.
#
# load_all(), not install_local(): nothing here runs on a crew worker.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_numbat_bulk_panel.R
