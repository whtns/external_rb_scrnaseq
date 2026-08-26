#!/bin/bash
#SBATCH --job-name=diag_filt_2clone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=00:50:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_filt_2clone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_filt_2clone_%j.err
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/diag_filtered_two_clone.R
