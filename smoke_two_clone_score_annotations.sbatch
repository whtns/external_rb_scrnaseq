#!/bin/bash
#SBATCH --job-name=smoke_two_clone_score_annot
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_two_clone_score_annot_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_two_clone_score_annot_%j.err

# Read-only against _targets_r431 (tar_read only, no tar_make) -- safe to run
# alongside nothing else; still check the queue before submitting.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_two_clone_score_annotations.R
