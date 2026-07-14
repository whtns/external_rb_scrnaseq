#!/bin/bash
#SBATCH --job-name=bigsur_diag
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/bigsur_diag_%j.log

# Usage: sbatch scripts/bigsur_cellranger.sbatch SRX10264523 [--correlations]

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

~/.conda/envs/ARMOR/bin/python scripts/bigsur_diagnostics.py "$@"
