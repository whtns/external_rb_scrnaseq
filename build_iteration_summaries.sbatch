#!/bin/bash
#SBATCH --job-name=iteration_summaries
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_iteration_summaries_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_iteration_summaries_%j.err

# Builds iteration_summaries: one PDF per eligible tumor with a column per numbat
# consensus round. 12 h rather than the 6 h the sample_summaries script uses --
# each of ~90 round-objects loads a joint_post over 100 MB, so this is the heavy
# one in the figure set.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# install_local, not load_all: these branches run on crew workers, which are fresh
# R sessions and never see the controller's load_all().
Rscript -e "devtools::install_local('/project2/cobrinik_1090/rpkgs/numbat_helpers', force = TRUE, upgrade = 'never')" \
  && Rscript -e "stopifnot(all(c('collate_iteration_summary','numbat_complete_rounds') %in% getNamespaceExports('numbatHelpers')))" \
  && Rscript run_targets_bg.R iteration_summaries iteration_summaries_report
