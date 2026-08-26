#!/bin/bash
#SBATCH --job-name=mosaic_interpret
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_interpret_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_interpret_%j.err

# Issue #40 interpretation add-on: the UNBIASED companion to the supplement.
# Reuses the completed production integration (output/mosaicmpi_rb/integrate_all)
# -- no cNMF re-run. Cheap: one broad ssGSEA + one R figure job.
#   1. build the broad HALLMARK-50 reference GMT (R)
#   2. score representative programs against it (mosaicmpi ssgsea, mosaic env)
#   3. assemble fig_mosaic_nmf_interpret.pdf + top-genes table + methods (R)
#   4. push the new deliverables to Drive (gdrive:rb_scrnaseq/mosaic_nmf/supplement)

set -euo pipefail
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

ROOT=output/mosaicmpi_rb
PKL="$ROOT/integrate_all/network_integration.pkl.gz"
HGMT=data/mosaic_hallmark.gmt
SAMPLES=(SRX10264523 SRX10264526 SRX11133592 SRX11133593 SRX11133594 SRX22868103 SRX10831287)
NW=$SLURM_CPUS_PER_TASK

R_ENV() {
  module purge; module load gcc/13.3.0; module load r/4.4.1
  export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}
MOSAIC_ENV() {
  module purge
  source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
  conda activate mosaic
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
         NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
}

echo "===== Stage 1: HALLMARK-50 GMT (R) ====="
R_ENV
Rscript src/mosaic/build_hallmark_gmt.R "$HGMT"

echo "===== Stage 2: broad ssGSEA (mosaic env) ====="
MOSAIC_ENV
mosaicmpi ssgsea -n "$PKL" -g "$HGMT" -o "$ROOT/ssgsea_hallmark" --cpus "$NW"

echo "===== Stage 3: interpretation figure + tables (R) ====="
R_ENV
Rscript src/mosaic/mosaic_interpret.R "$ROOT/ssgsea_hallmark" "$ROOT/usage" "$ROOT/integrate_all" "${SAMPLES[@]}"

echo "===== Stage 4: push to Drive ====="
module load rclone
DEST=gdrive:rb_scrnaseq/mosaic_nmf/supplement
rclone copyto results/fig_mosaic_nmf_interpret.pdf        "$DEST/fig_mosaic_nmf_interpret.pdf"
rclone copyto results/table_mosaic_nmf_community_genes.csv "$DEST/table_mosaic_nmf_community_genes.csv"
rclone copyto docs/mosaic_nmf_interpret.md                "$DEST/mosaic_nmf_interpret.md"
rclone copyto "$ROOT/ssgsea_hallmark/representative_program_nes.pdf" "$DEST/ssgsea_hallmark_nes.pdf"

echo "INTERPRET_OK: figure=$(ls -la results/fig_mosaic_nmf_interpret.pdf 2>/dev/null | awk '{print $5}') bytes"

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"mosaic_interpret (#40 add-on) $SLURM_JOB_ID finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
