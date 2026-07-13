#!/bin/bash
#SBATCH --job-name=eval_hyp_phase
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/eval_hypoxia_cluster_phase_%j.log

set -euo pipefail
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e "devtools::load_all('/project2/cobrinik_1090/rpkgs/numbat_helpers'); source('src/eval_hypoxia_cluster_phase.R')"
