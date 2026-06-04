#!/bin/bash
#SBATCH --job-name=fig_s03a_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH -o build_fig_s03a_load_all.log

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e "devtools::load_all('/project2/cobrinik_1090/rpkgs/numbat_helpers')" && \
  Rscript run_targets_bg.R fig_s03a_unfiltered_plots
