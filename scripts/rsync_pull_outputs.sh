#!/usr/bin/env bash
# Run this script ON HPC to pull outputs from the source machine — avoids 2FA
# since HPC initiates all connections back to SOURCE_HOST (no inbound auth needed).
#
# Usage (on HPC):
#   [SOURCE=user@host] [COPY_FASTQ=true] bash scripts/rsync_pull_outputs.sh
#
# Requires SSH key access from HPC to SOURCE_HOST.

set -euo pipefail

SOURCE_HOST=${SOURCE_HOST:-skevin@cobrinik-1.saban-chla.usc.edu}
SOURCE_BASE=/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj
DEST=/project2/cobrinik_1090/external_rb_scrnaseq_proj
COPY_FASTQ=${COPY_FASTQ:-false}

RSYNC="rsync -rltvh --append-verify --progress"

SAMPLES=$(ssh "$SOURCE_HOST" "ls '${SOURCE_BASE}/data/FASTQ' | grep '^SRR'")

# ── cellranger filtered matrices ──────────────────────────────────────────────
for s in $SAMPLES; do
    src="${SOURCE_HOST}:${SOURCE_BASE}/output/cellranger/${s}/outs/filtered_feature_bc_matrix"
    dst="${DEST}/output/cellranger/${s}/outs/filtered_feature_bc_matrix"
    if ssh "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/output/cellranger/${s}/outs/filtered_feature_bc_matrix' ]]"; then
        mkdir -p "$dst"
        $RSYNC "${src}/" "${dst}/"
    fi
done

# ── numbat filtered results ───────────────────────────────────────────────────
for s in $SAMPLES; do
    src="${SOURCE_HOST}:${SOURCE_BASE}/output/numbat_sridhar_filtered/${s}"
    dst="${DEST}/output/numbat_sridhar_filtered/${s}"
    if ssh "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/output/numbat_sridhar_filtered/${s}' ]]"; then
        mkdir -p "$dst"
        $RSYNC "${src}/" "${dst}/"
    fi
done

# ── per-sample seurat objects (filtered + regressed) ─────────────────────────
mkdir -p "${DEST}/output/seurat"
for s in $SAMPLES; do
    for pattern in "${s}_filtered_seu.rds" "${s}_regressed_seu.rds"; do
        src="${SOURCE_HOST}:${SOURCE_BASE}/output/seurat/${pattern}"
        if ssh "$SOURCE_HOST" "[[ -f '${SOURCE_BASE}/output/seurat/${pattern}' ]]"; then
            $RSYNC "$src" "${DEST}/output/seurat/"
        fi
    done
done

# ── fastq files (opt-in: COPY_FASTQ=true) ────────────────────────────────────
if [[ "$COPY_FASTQ" == "true" ]]; then
    for s in $SAMPLES; do
        src="${SOURCE_HOST}:${SOURCE_BASE}/data/FASTQ/${s}"
        dst="${DEST}/data/FASTQ/${s}"
        if ssh "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/data/FASTQ/${s}' ]]"; then
            mkdir -p "$dst"
            $RSYNC "${src}/" "${dst}/"
        fi
    done
fi

# ── targets store ─────────────────────────────────────────────────────────────
mkdir -p "${DEST}/_targets_r431"
$RSYNC "${SOURCE_HOST}:${SOURCE_BASE}/_targets_r431/" "${DEST}/_targets_r431/"
