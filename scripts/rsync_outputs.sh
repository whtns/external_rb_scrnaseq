#!/usr/bin/env bash
# Copy critical snakemake and targets outputs for all samples in data/FASTQ.
# Usage: [DEST=user@host:/path] [COPY_FASTQ=true] bash scripts/rsync_outputs.sh

set -euo pipefail

PROJ=/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj
DEST=${DEST:-stachele@hpc-transfer1.usc.edu:/project2/cobrinik_1090/external_rb_scrnaseq_proj}
COPY_FASTQ=${COPY_FASTQ:-false}

SAMPLES=$(ls "${PROJ}/data/FASTQ" | grep '^SRR')

# -r recursive  -l copy symlinks  -t preserve times  -v verbose  -h human sizes
# --append-verify resumes interrupted transfers and re-checks checksum on completion
RSYNC="rsync -rltvh --append-verify --progress"

# Extract host and base path from DEST for pre-creating directories
_DEST_HOST="${DEST%%:*}"
_DEST_BASE="${DEST#*:}"

_ssh_mkdir() {
    # $1 = subpath relative to _DEST_BASE
    if [[ "$DEST" == *":"* ]]; then
        ssh "$_DEST_HOST" "mkdir -p '${_DEST_BASE}/$1'"
    else
        mkdir -p "${_DEST_BASE}/$1"
    fi
}

# ── cellranger filtered matrices (skip BAMs) ──────────────────────────────────
for s in $SAMPLES; do
    src="${PROJ}/output/cellranger/${s}/outs/filtered_feature_bc_matrix"
    if [[ -d "$src" ]]; then
        _ssh_mkdir "output/cellranger/${s}/outs/filtered_feature_bc_matrix"
        $RSYNC "$src/" "${DEST}/output/cellranger/${s}/outs/filtered_feature_bc_matrix/"
    fi
done

# ── numbat filtered results ───────────────────────────────────────────────────
for s in $SAMPLES; do
    src="${PROJ}/output/numbat_sridhar_filtered/${s}"
    if [[ -d "$src" ]]; then
        _ssh_mkdir "output/numbat_sridhar_filtered/${s}"
        $RSYNC "$src/" "${DEST}/output/numbat_sridhar_filtered/${s}/"
    fi
done

# ── per-sample seurat objects (filtered + regressed) ─────────────────────────
_ssh_mkdir "output/seurat"
for s in $SAMPLES; do
    for pattern in "${s}_filtered_seu.rds" "${s}_regressed_seu.rds"; do
        src="${PROJ}/output/seurat/${pattern}"
        [[ -f "$src" ]] && $RSYNC "$src" "${DEST}/output/seurat/"
    done
done

# ── fastq files (opt-in: COPY_FASTQ=true) ────────────────────────────────────
if [[ "$COPY_FASTQ" == "true" ]]; then
    for s in $SAMPLES; do
        src="${PROJ}/data/FASTQ/${s}"
        if [[ -d "$src" ]]; then
            _ssh_mkdir "data/FASTQ/${s}"
            $RSYNC "$src/" "${DEST}/data/FASTQ/${s}/"
        fi
    done
fi

# ── targets store ─────────────────────────────────────────────────────────────
_ssh_mkdir "_targets_r431"
$RSYNC "${PROJ}/_targets_r431/" "${DEST}/_targets_r431/"
