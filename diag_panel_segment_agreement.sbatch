#!/bin/bash
#SBATCH --job-name=diag_seg_agree
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_seg_agree_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_seg_agree_%j.err

# Does the #42 fix actually make the heatmap and the bulk panel show the same
# segments? Reads the numbat RDS files only -- no targets store, safe to run
# alongside the build.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript src/diag_panel_segment_agreement.R
