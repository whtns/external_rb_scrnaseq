#!/bin/bash
#SBATCH --job-name=rclone_hypoxia_rebuilt
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rclone_hypoxia_rebuilt_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rclone_hypoxia_rebuilt_%j.err

# Build the hypoxia_rebuilt_gdrive target: rclone the rebuilt scope-b hypoxia
# deliverables (per-sample hypoxia summaries, low-hypoxia clone/segment trees,
# low-hypoxia numbat + hypoxia-gene heatmaps, annotated marker collages) to
# gdrive:rb_scrnaseq/hypoxia_rebuilt.
#
# The target is tagged deployment = "main", so it runs in THIS process rather
# than on a crew worker -- rclone (loaded below) must therefore be on PATH here.
# All upstream inputs are already cached from the completed rebuild, so tar_make
# only runs the rclone copy.

module load r/4.4.1 rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript run_targets_bg.R hypoxia_rebuilt_gdrive
