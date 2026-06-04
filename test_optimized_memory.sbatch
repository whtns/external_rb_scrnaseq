#!/bin/bash
#SBATCH --job-name=test_memory_opt
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --output=test_memory_opt_%j.log
#SBATCH --error=test_memory_opt_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Test with a smaller target to verify everything works
Rscript run_targets_bg.R sample_summaries
