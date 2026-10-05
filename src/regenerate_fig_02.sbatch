#!/bin/bash
#SBATCH --job-name=regen-fig02
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/regen_fig02_%j.log

module load r/4.4.1 poppler/23.04.0
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/regenerate_fig_02.R
