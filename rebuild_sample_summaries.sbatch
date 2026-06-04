#!/bin/bash
#SBATCH --job-name=sample_summaries
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --output=rebuild_sample_summaries_%j.log
#SBATCH --error=rebuild_sample_summaries_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript run_targets_bg.R sample_summaries
