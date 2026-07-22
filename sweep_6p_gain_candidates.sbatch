#!/bin/bash
#SBATCH --job-name=sweep_6p_gain
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/sweep_6p_gain_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/sweep_6p_gain_%j.err

# Tier 1 sweep: 6p gains that fell below min_LLR = 5, SRX samples only (39 dirs).
# Read-only over output/numbat_sridhar; touches no targets store, no rerun.
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/sweep_6p_gain_candidates.R
