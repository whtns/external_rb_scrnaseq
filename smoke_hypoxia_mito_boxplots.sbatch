#!/bin/bash
#SBATCH --job-name=smoke_hyp_mito_box
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_hyp_mito_box_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_hyp_mito_box_%j.err

# Smoke test only -- three samples, no targets store touched, writes
# results/hypoxia_cluster_split/SMOKE_hypoxia_mito_score_boxplots.pdf.
# load_all(), not install_local(): nothing here runs on a crew worker.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_hypoxia_mito_boxplots.R
