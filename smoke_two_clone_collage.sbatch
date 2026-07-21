#!/bin/bash
#SBATCH --job-name=smoke_two_clone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_two_clone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_two_clone_%j.err

# Smoke test the new plot_scna_two_clone_res_collages() on one paper_retained
# sample before committing to the full sweep build.
module load r/4.4.1
module load poppler/23.04.0
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_two_clone_collage.R
