#!/bin/bash
#SBATCH --job-name=fig_s04
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --output=build_fig_s04_%j.log
#SBATCH --error=build_fig_s04_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript run_targets_bg.R fig_s04
