#!/bin/bash
#SBATCH --job-name=diag_convergence
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_convergence_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_convergence_%j.err

# How many numbat iterations were actually needed? Measured from the per-round
# segs_consensus_<k>.tsv files already on disk, all 71 samples, no rerun.
# Read-only; touches no targets store.
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/diag_consensus_convergence.R
