#!/bin/bash
#SBATCH --job-name=smoke_collage_panels
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_collage_panels_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_collage_panels_%j.err

# Smoke test for issues #34/#35/#37/#38 (bar colours, diagonal heatmap labels,
# segment-tree panel, clustree panel). Uses devtools::load_all so the smoke runs
# the working-tree code without an install; a single-process test, no crew.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_collage_panels.R
