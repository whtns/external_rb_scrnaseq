#!/bin/bash
#SBATCH --job-name=test_fig_s03a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH -o test_fig_s03a_%j.log

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS="/project2/cobrinik_1090/external_rb_scrnaseq_proj/.R_libs:${R_LIBS}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Run the target with run_targets_bg.R (which now has devtools::load_all())
Rscript run_targets_bg.R fig_s03a_unfiltered_plots
