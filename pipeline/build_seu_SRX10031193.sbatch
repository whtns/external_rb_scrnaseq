#!/bin/bash
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=40:00:00
#SBATCH --job-name=seu-SRX10031193
#SBATCH --output=logs/slurm/build_seu_SRX10031193-%j.out
#SBATCH --error=logs/slurm/build_seu_SRX10031193-%j.err

set -euo pipefail

module load apptainer

# snakemake lives in the ARMOR conda env
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate ARMOR
export PATH="$HOME/.local/bin:$PATH"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj/pipeline

TARGET=/project2/cobrinik_1090/external_rb_scrnaseq_proj/output/seurat/SRX10031193_seu.rds

# Builds the cellranger -> seurat chain only (no pileup/phase/numbat).
snakemake \
    --cores 32 \
    --rerun-incomplete \
    --printshellcmds \
    "$TARGET"
