#!/bin/bash
#SBATCH --job-name=diag_hla_6p
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_hla_6p_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_hla_6p_%j.err

# Issue #31, step 0: quantify HLA/MHC exposure of every called chr6p segment from
# the numbat output we already have, before committing to any masked rerun.
# Read-only over output/numbat_sridhar; touches no targets store.
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/diag_hla_6p_contribution.R
