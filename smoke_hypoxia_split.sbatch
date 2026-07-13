#!/bin/bash
#SBATCH --job-name=smoke_hypoxia_split
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/smoke_hypoxia_split_%j.log

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/smoke_hypoxia_split.R
