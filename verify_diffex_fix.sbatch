#!/bin/bash
#SBATCH --job-name=verify_diffex_fix
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/verify_diffex_fix_%j.log

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/verify_diffex_cluster_fix.R
