#!/bin/bash
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --job-name=numbat-SRX10031193
#SBATCH --output=logs/slurm/build_numbat_SRX10031193-%j.out
#SBATCH --error=logs/slurm/build_numbat_SRX10031193-%j.err

set -euo pipefail

module load apptainer

# snakemake lives in the ARMOR conda env
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate ARMOR
export PATH="$HOME/.local/bin:$PATH"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj/pipeline

TARGET=/project2/cobrinik_1090/external_rb_scrnaseq_proj/output/numbat_sridhar/SRX10031193_numbat.rds

# Builds pileup -> phase_and_count -> numbat_sridhar -> numbat_sridhar_rds.
# num_clone=5 (matches sibling collin samples SRX10031191/192/194).
snakemake \
    --cores 16 \
    --rerun-incomplete \
    --printshellcmds \
    "$TARGET"
