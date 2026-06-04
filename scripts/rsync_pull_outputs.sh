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

# Reuse one SSH connection for all ssh/rsync calls — one password prompt only
SSH_SOCKET="/tmp/rsync_pull_$$.sock"
ssh -fNM -o ControlMaster=yes -o ControlPath="$SSH_SOCKET" \
    -o ControlPersist=10m "$SOURCE_HOST"
trap 'ssh -O exit -o ControlPath="$SSH_SOCKET" "$SOURCE_HOST" 2>/dev/null; rm -f "$SSH_SOCKET"' EXIT

SSH=(ssh -o "ControlPath=$SSH_SOCKET")
RSYNC=(rsync -rltvh --append-verify --progress -e "ssh -o ControlPath=$SSH_SOCKET")

SAMPLES=$("${SSH[@]}" "$SOURCE_HOST" "ls '${SOURCE_BASE}/data/FASTQ' | grep '^SRR'")

# ── cellranger BAM + index (large files — resumable via --append-verify) ──────
# (already copied and renamed to SRX IDs on HPC — uncomment to re-pull)
# for s in $SAMPLES; do
#     outs_src="${SOURCE_HOST}:${SOURCE_BASE}/output/cellranger/${s}/outs"
#     outs_dst="${DEST}/output/cellranger/${s}/outs"
#     for f in possorted_genome_bam.bam possorted_genome_bam.bam.bai; do
#         if "${SSH[@]}" "$SOURCE_HOST" "[[ -f '${SOURCE_BASE}/output/cellranger/${s}/outs/${f}' ]]"; then
#             mkdir -p "$outs_dst"
#             "${RSYNC[@]}" "${outs_src}/${f}" "${outs_dst}/"
#         fi
#     done
# done

# ── cellranger filtered matrices ──────────────────────────────────────────────
# (already copied and renamed to SRX IDs on HPC — uncomment to re-pull)
# for s in $SAMPLES; do
#     src="${SOURCE_HOST}:${SOURCE_BASE}/output/cellranger/${s}/outs/filtered_feature_bc_matrix"
#     dst="${DEST}/output/cellranger/${s}/outs/filtered_feature_bc_matrix"
#     if "${SSH[@]}" "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/output/cellranger/${s}/outs/filtered_feature_bc_matrix' ]]"; then
#         mkdir -p "$dst"
#         "${RSYNC[@]}" "${src}/" "${dst}/"
#     fi
# done

# ── per-sample seurat objects (filtered + regressed) ─────────────────────────
# (already copied and renamed to SRX IDs on HPC — uncomment to re-pull)
# mkdir -p "${DEST}/output/seurat"
# for s in $SAMPLES; do
#     for pattern in "${s}_filtered_seu.rds" "${s}_regressed_seu.rds"; do
#         src="${SOURCE_HOST}:${SOURCE_BASE}/output/seurat/${pattern}"
#         if "${SSH[@]}" "$SOURCE_HOST" "[[ -f '${SOURCE_BASE}/output/seurat/${pattern}' ]]"; then
#             "${RSYNC[@]}" "$src" "${DEST}/output/seurat/"
#         fi
#     done
# done

# ── per-sample base seurat objects (_seu.rds) ─────────────────────────────────
# Source machine names files by SRR ID; HPC uses SRX IDs for renamed samples.
# Reads SRX->SRR mapping from metadata at runtime and renames on copy.
mkdir -p "${DEST}/output/seurat"

# Build associative array: SRX -> SRR (where Run=SRR and Experiment=SRX)
declare -A SRX_TO_SRR
while IFS=$'\t' read -r run _ _ _ _ _ exp rest; do
    if [[ "$run" == SRR* && "$exp" == SRX* ]]; then
        SRX_TO_SRR["$exp"]="$run"
    fi
done < <(tail -n +2 "${DEST}/data/metadata.tsv")

for s in $(ls "${DEST}/output/cellranger/"); do
    [[ -f "${DEST}/output/seurat/${s}_seu.rds" ]] && continue
    # SRX samples that were renamed from an SRR on the source machine
    src_name="$s"
    if [[ "$s" == SRX* && -n "${SRX_TO_SRR[$s]+x}" ]]; then
        src_name="${SRX_TO_SRR[$s]}"
    fi
    if "${SSH[@]}" "$SOURCE_HOST" "[[ -f '${SOURCE_BASE}/output/seurat/${src_name}_seu.rds' ]]"; then
        "${RSYNC[@]}" \
            "${SOURCE_HOST}:${SOURCE_BASE}/output/seurat/${src_name}_seu.rds" \
            "${DEST}/output/seurat/${s}_seu.rds"
    fi
done

# ── numbat filtered results ───────────────────────────────────────────────────
for s in $SAMPLES; do
    src="${SOURCE_HOST}:${SOURCE_BASE}/output/numbat_sridhar_filtered/${s}"
    dst="${DEST}/output/numbat_sridhar_filtered/${s}"
    if "${SSH[@]}" "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/output/numbat_sridhar_filtered/${s}' ]]"; then
        mkdir -p "$dst"
        "${RSYNC[@]}" "${src}/" "${dst}/"
    fi
done

# ── numbat_sridhar results directory (contains done.txt + all numbat outputs) ─
for s in $SAMPLES; do
    src="${SOURCE_HOST}:${SOURCE_BASE}/output/numbat_sridhar/${s}"
    dst="${DEST}/output/numbat_sridhar/${s}"
    if "${SSH[@]}" "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/output/numbat_sridhar/${s}' ]]"; then
        mkdir -p "$dst"
        "${RSYNC[@]}" "${src}/" "${dst}/"
    fi
done

# ── numbat_sridhar packed RDS ─────────────────────────────────────────────────
mkdir -p "${DEST}/output/numbat_sridhar"
for s in $SAMPLES; do
    src="${SOURCE_HOST}:${SOURCE_BASE}/output/numbat_sridhar/${s}_numbat.rds"
    if "${SSH[@]}" "$SOURCE_HOST" "[[ -f '${SOURCE_BASE}/output/numbat_sridhar/${s}_numbat.rds' ]]"; then
        "${RSYNC[@]}" "$src" "${DEST}/output/numbat_sridhar/"
    fi
done

# ── fastq files (opt-in: COPY_FASTQ=true) ────────────────────────────────────
if [[ "$COPY_FASTQ" == "true" ]]; then
    for s in $SAMPLES; do
        src="${SOURCE_HOST}:${SOURCE_BASE}/data/FASTQ/${s}"
        dst="${DEST}/data/FASTQ/${s}"
        if "${SSH[@]}" "$SOURCE_HOST" "[[ -d '${SOURCE_BASE}/data/FASTQ/${s}' ]]"; then
            mkdir -p "$dst"
            "${RSYNC[@]}" "${src}/" "${dst}/"
        fi
    done
fi

# ── targets store ─────────────────────────────────────────────────────────────
mkdir -p "${DEST}/_targets_r431"
"${RSYNC[@]}" "${SOURCE_HOST}:${SOURCE_BASE}/_targets_r431/" "${DEST}/_targets_r431/"
