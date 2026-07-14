#!/bin/bash
#SBATCH --job-name=bigsur_all
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=1-39%10
#SBATCH --output=logs/bigsur_all_%A_%a.log

# BigSur feature selection + correlations + diagnostics for every sample with a
# CellRanger filtered matrix (see scripts/bigsur_samples.txt).
#
# Baseline (no gene filtering):
#   sbatch scripts/bigsur_all_samples.sbatch
#
# Dropping the haemoglobin / crystallin / immunoglobulin carry-over families:
#   sbatch --export=ALL,BIGSUR_OUTDIR=results/bigsur_nocontam,BIGSUR_EXCLUDE=1 \
#          scripts/bigsur_all_samples.sbatch

set -euo pipefail

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/bigsur_samples.txt)
OUTDIR=${BIGSUR_OUTDIR:-results/bigsur}
EXCLUDE=${BIGSUR_EXCLUDE:-}
echo "=== [${SLURM_ARRAY_TASK_ID}] ${SAMPLE} -> ${OUTDIR} ==="

export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
PY=~/.conda/envs/ARMOR/bin/python

# Bare --exclude-contaminants drops all families; both scripts must agree on it
# or the diagnostic panels won't match the correlation matrix they plot.
EXCL_ARG=()
[[ -n "$EXCLUDE" ]] && EXCL_ARG=(--exclude-contaminants)

$PY scripts/bigsur_cellranger.py "$SAMPLE" --correlations --n-corr-genes 2000 \
    --outdir "$OUTDIR" "${EXCL_ARG[@]}"
$PY scripts/bigsur_diagnostics.py "$SAMPLE" --outroot "$OUTDIR" "${EXCL_ARG[@]}"

echo "=== DONE ${SAMPLE} ==="
