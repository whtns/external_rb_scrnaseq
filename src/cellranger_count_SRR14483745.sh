#!/bin/bash
#SBATCH --job-name=cr_SRR14483745
#SBATCH --output=logs/cellranger_SRR14483745_%j.out
#SBATCH --error=logs/cellranger_SRR14483745_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

SAMPLE=SRR14483745
FASTQ_DIR=/project2/cobrinik_1090/external_rb_scrnaseq_proj/data/FASTQ/${SAMPLE}
TRANSCRIPTOME=/project2/cobrinik_1090/Homo_sapiens/cellranger/refdata-gex-GRCh38-2024-A
OUTDIR=/project2/cobrinik_1090/external_rb_scrnaseq_proj/output/cellranger

mkdir -p logs "${OUTDIR}"

module load cellranger/8.0.1

cd "${OUTDIR}"

cellranger count \
    --id="${SAMPLE}" \
    --sample="${SAMPLE}" \
    --fastqs="${FASTQ_DIR}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --chemistry=SC3Pv2 \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=120 \
    --create-bam=true
