#!/bin/bash
#SBATCH --job-name=diag_diffex_bug
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_diffex_bug_%j.log

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/diag_diffex_cluster_bug.R
