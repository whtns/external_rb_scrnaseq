#!/bin/bash
#SBATCH --job-name=smoke_bar_signif
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_bar_signif_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_bar_signif_%j.err

# Issue #32: per-cluster significance marks on the collage stacked-bar panel.
# Verify the marks, the caption, and the sidecar enrichment CSV on one
# sample/SCNA/resolution before rebuilding the whole two-clone sweep.
# Uses devtools::load_all, so this tests the WORKING TREE, not the installed pkg.
module load r/4.4.1
module load poppler/23.04.0
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_bar_signif.R
