#!/bin/bash
#SBATCH --job-name=vfy_rds
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/vfy_rds_%j.log

# Targeted five-locus test for canonical RB SCNAs over the 39 SRX samples.
# Reads bulk_clones_final.tsv.gz only -- no numbat rerun, no targets store, so
# this is safe to run alongside anything else.

set -euo pipefail

module purge
# r/4.4.1 lives under gcc/13.3.0 in the Lmod hierarchy -- after a purge it is
# not visible until the compiler is loaded (module spider r/4.4.1).
module load gcc/13.3.0
module load r/4.4.1

export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/verify_rebuilt_rds.R
