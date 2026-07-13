#!/bin/bash
#SBATCH --job-name=diag_hypoxia_low_scna
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:40:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_hypoxia_low_scna_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_hypoxia_low_scna_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/diag_hypoxia_low_scna.R
