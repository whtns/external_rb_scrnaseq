#!/bin/bash
#SBATCH --job-name=dl_prjna728722
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/src/logs/download_%j.out
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/src/logs/download_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --partition=main

set -euo pipefail

mkdir -p logs

module load sratoolkit/3.2.1

# Redirect SRA cache away from home directory (limited quota)
export NCBI_SETTINGS=/project2/cobrinik_1090/external_rb_scrnaseq_proj/ncbi/user-settings.mkfg
mkdir -p /project2/cobrinik_1090/external_rb_scrnaseq_proj/ncbi/sra_cache
vdb-config -s /repository/user/main/public/root=/project2/cobrinik_1090/external_rb_scrnaseq_proj/ncbi/sra_cache 2>/dev/null || true

OUTDIR="/project2/cobrinik_1090/external_rb_scrnaseq_proj/fastq"
SCRIPT_DIR="/project2/cobrinik_1090/external_rb_scrnaseq_proj/src"

python3 "${SCRIPT_DIR}/download_prjna728722.py" \
    --outdir "$OUTDIR" \
    --threads 32 \
    --batch-size 50
