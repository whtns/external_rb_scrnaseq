#!/bin/bash
#SBATCH --job-name=r_pak_restore
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_restore_sbatch_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_restore_sbatch_%j.log

set -euo pipefail

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load R module
module load r/4.4.1

# Set up LD_LIBRARY_PATH for libRlapack.so/libRblas.so/libhdf5 symlinks
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Source-only mode; all tarballs pre-cached from login/OnDemand node — no internet needed
export PKGCACHE_PLATFORMS="x86_64-pc-linux-gnu"
export R_PKG_CACHE_DIR="$HOME/.cache/R/pkgcache"
# Prevent metadata refresh so pak uses cached PACKAGES files without internet
export PKGCACHE_UPDATE_AFTER=Inf

# Run pak restore — installs to ~/R (not .R_libs) because lib is explicitly set
Rscript /project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_restore_full.R

echo "Job finished: $(date)"
