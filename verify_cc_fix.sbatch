#!/bin/bash
#SBATCH --job-name=verify_cc_fix
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=00:40:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/verify_cc_fix_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/verify_cc_fix_%j.err
module load gcc/13.3.0; module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/verify_cc_fix.R
