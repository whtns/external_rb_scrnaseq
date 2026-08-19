#!/bin/bash
#SBATCH --job-name=diag_hyp_mito_cor
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_hyp_mito_cor_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_hyp_mito_cor_%j.err

# Read-only diagnostic: per-cell correlation between the two module scores.
# Touches no targets store and writes one CSV.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/diag_hypoxia_mito_correlation.R
