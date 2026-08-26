#!/bin/bash
#SBATCH --job-name=rb_segment_trees
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rb_segment_trees_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rb_segment_trees_%j.err

# Segment-labeled clone trees for the samples in
# results/rb_scna_gain_by_round_selection.csv. Read-only with respect to
# output/numbat_sridhar/; writes only under results/. No tar_make, so the
# one-job rule does not apply.
set -euo pipefail

module load gcc/13.3.0 r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/plot_gained_segment_trees.R
