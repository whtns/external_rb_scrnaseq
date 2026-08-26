#!/bin/bash
#SBATCH --job-name=verify_maj
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=logs/verify_majority_%j.log

# Verify every *_numbat.rds holds the round results/numbat_selected_round.csv
# now specifies, and carries a usable phylogeny. Reads 36 objects, several of
# them >150 MB, so it needs a real allocation rather than the login node.

set -euo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
cd "$PROJ"

module purge
module load gcc/13.3.0 r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

Rscript src/verify_majority_rds.R
