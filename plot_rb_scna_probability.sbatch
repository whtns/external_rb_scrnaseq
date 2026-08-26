#!/bin/bash
#SBATCH --job-name=rb_scna_prob
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rb_scna_prob_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rb_scna_prob_%j.err

# Per-cell SCNA posterior plots at the canonical RB events for the samples in
# results/rb_scna_gain_by_round_selection.csv. Read-only with respect to
# output/numbat_sridhar/; writes only under results/. Does NOT run tar_make, so
# the one-job rule does not apply.
set -euo pipefail

module load gcc/13.3.0 r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# plot_rb_scna_probability() is new; make sure the installed package has it.
Rscript -e "devtools::install_local('/project2/cobrinik_1090/rpkgs/numbat_helpers', force = TRUE, upgrade = 'never')"

Rscript src/plot_rb_scna_probability.R
