#!/bin/bash
#SBATCH --job-name=unfiltered_seus_build
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --output=rebuild_unfiltered_seus_%j.log
#SBATCH --error=rebuild_unfiltered_seus_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript run_targets_bg.R unfiltered_seus
