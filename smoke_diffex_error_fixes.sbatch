#!/bin/bash
#SBATCH --job-name=smoke_diffex_fixes
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_diffex_fixes_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_diffex_fixes_%j.err

module load r/4.4.1
module load poppler/23.04.0
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_diffex_error_fixes.R
