#!/bin/bash
#SBATCH --job-name=dl_prjna728722
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/src/logs/download_%A_%a.out
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/src/logs/download_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=main
#SBATCH --array=0-56

set -euo pipefail

# All 57 target accessions (index matches $SLURM_ARRAY_TASK_ID)
ACCESSIONS=(
    SRR14483745 SRR14483746 SRR14483747 SRR14483748 SRR14483749
    SRR14483750 SRR14483751 SRR14483752 SRR14483753 SRR14483754
    SRR14483755 SRR14483756 SRR14483757 SRR14483758 SRR14483759
    SRR14483760 SRR14483761 SRR14483762 SRR14483763 SRR14483764
    SRR14483765 SRR14483766 SRR14483767 SRR14483768 SRR14483769
    SRR14483770 SRR14483771 SRR14483772 SRR14483773 SRR14483774
    SRR14483775 SRR14483776
    SRR14483781 SRR14483782 SRR14483783 SRR14483784
    SRR14483937 SRR14483938 SRR14483939
    SRR14483958 SRR14483959 SRR14483960 SRR14483961 SRR14483962
    SRR14483963 SRR14483964 SRR14483965 SRR14483966 SRR14483967
    SRR14483968 SRR14483969 SRR14483970 SRR14483971 SRR14483972
    SRR14483973 SRR14483974 SRR14483975
)

ACC="${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}"
OUTDIR="/project2/cobrinik_1090/external_rb_scrnaseq_proj/fastq"
TMPDIR="/project2/cobrinik_1090/external_rb_scrnaseq_proj/fastq/tmp/${ACC}"
SRA_CACHE="/project2/cobrinik_1090/external_rb_scrnaseq_proj/ncbi/sra_cache"
NCBI_SETTINGS="/project2/cobrinik_1090/external_rb_scrnaseq_proj/ncbi/user-settings.mkfg"

mkdir -p "$OUTDIR" "$TMPDIR" "$SRA_CACHE"

module load sratoolkit/3.2.1

export NCBI_SETTINGS
if [ ! -f "$NCBI_SETTINGS" ]; then
    vdb-config -s /repository/user/main/public/root="$SRA_CACHE"
fi

# Skip if all expected files already exist and are non-empty
existing=$(find "$OUTDIR" -name "${ACC}_*.fastq.gz" -size +0c 2>/dev/null | wc -l)
if [ "$existing" -gt 0 ]; then
    echo "[$SLURM_ARRAY_TASK_ID] $ACC already has $existing file(s), skipping"
    exit 0
fi

echo "[$SLURM_ARRAY_TASK_ID] fasterq-dump $ACC ..."

fasterq-dump "$ACC" \
    --outdir "$OUTDIR" \
    --temp "$TMPDIR" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --split-files \
    --include-technical \
    --progress

# Compress in parallel
pigz -p "$SLURM_CPUS_PER_TASK" "$OUTDIR/${ACC}"_*.fastq

# Clean up per-accession temp dir
rm -rf "$TMPDIR"

echo "[$SLURM_ARRAY_TASK_ID] Done: $ACC"
find "$OUTDIR" -name "${ACC}_*.fastq.gz" -size +0c
