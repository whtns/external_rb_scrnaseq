#!/bin/bash
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --job-name=install_nbh
#SBATCH --output=logs/install_nbh_%j.log
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript -e "devtools::install_local('/project2/cobrinik_1090/rpkgs/numbat_helpers', force = TRUE, upgrade = 'never')" && echo INSTALL_OK
