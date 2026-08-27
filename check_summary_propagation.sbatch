#!/bin/bash
#SBATCH --job-name=chk_prop
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/chk_prop_%j.log

# Read-only dependency-propagation check for sample_summaries. Builds nothing.
set -euo pipefail
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
module purge
module load gcc/13.3.0 r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"
Rscript src/check_summary_propagation.R
