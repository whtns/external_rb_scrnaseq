#!/bin/bash
#SBATCH --job-name=mosaic_supplement
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=24:00:00
#SBATCH --mem=160G
#SBATCH --cpus-per-task=32
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_supplement_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_supplement_%j.err

# Issue #40 PRODUCTION: full mosaicMPI NMF supplement over the 7 paper-retained
# tumors (A->G). Validated end-to-end at small scale by rebuild_mosaic_smoke.sbatch.
# Deliverables: results/fig_mosaic_nmf_supplement.pdf, results/table_mosaic_nmf_programs.csv,
# docs/mosaic_nmf_supplement.md.
#
# cNMF is the cost: k=2..15 x 100 iter x 7 samples. factorize is fanned out across
# $NW single-threaded workers (thread pinning below is essential -- multi-threaded
# BLAS x many procs oversubscribes the node and is ~16x slower). ~1 min/task /
# 32 workers -> a few hours.

set -euo pipefail
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

SAMPLES=(SRX10264523 SRX10264526 SRX11133592 SRX11133593 SRX11133594 SRX22868103 SRX10831287)
ROOT=output/mosaicmpi_rb
GMT=data/mosaic_reference_programs.gmt
KFIRST=2; KLAST=15
NITER=100
SEED=42
NGENES=2000
NW=$SLURM_CPUS_PER_TASK
mkdir -p logs "$ROOT/inputs"

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

echo "===== Stage A + counts export (R) ====="
R_ENV
Rscript src/mosaic/build_reference_gmt.R "$GMT"
for s in "${SAMPLES[@]}"; do Rscript src/mosaic/export_counts.R "$s" "$ROOT/$s"; done

echo "===== Stage B-C factorize/postprocess (mosaic env) ====="
MOSAIC_ENV
for s in "${SAMPLES[@]}"; do
  echo "--- $s ---"
  python src/mosaic/build_h5ad.py "$ROOT/$s" "$ROOT/inputs/${s}.raw.h5ad"
  mosaicmpi check-h5ad     -i "$ROOT/inputs/${s}.raw.h5ad" -o "$ROOT/inputs/${s}.h5ad"
  mosaicmpi model-odg      -n "$s" -o "$ROOT" -i "$ROOT/inputs/${s}.h5ad"
  python src/mosaic/set_parameters.py "$ROOT" "$s" "$NGENES" $KFIRST $KLAST "$NITER" "$SEED"
  for i in $(seq 0 $((NW-1))); do
    mosaicmpi factorize -n "$s" -o "$ROOT" --worker_index "$i" --total_workers "$NW" &
  done
  wait
  # local_density_threshold 0.5: filter outlier replicate spectra for cleaner
  # consensus programs (2.0 = no filtering). Revisit from the k-selection diagnostics.
  mosaicmpi postprocess    -n "$s" -o "$ROOT" --cpus "$NW" --local_density_threshold 0.5
done

echo "===== Stage D create-config + integrate ====="
CFG_ARGS=(); for s in "${SAMPLES[@]}"; do CFG_ARGS+=(-i "$ROOT/$s/${s}.h5ad"); done
mosaicmpi create-config "${CFG_ARGS[@]}" -o "$ROOT/config.toml"
mosaicmpi integrate -o "$ROOT/integrate_all" -c "$ROOT/config.toml" --cpus "$NW"
PKL=$(ls "$ROOT/integrate_all"/network_integration.pkl* 2>/dev/null | head -1)
echo "integration pkl: $PKL"

echo "===== Stage E ssgsea ====="
mosaicmpi ssgsea -n "$PKL" -g "$GMT" -o "$ROOT/ssgsea" --cpus "$NW"

echo "===== Stage F community usage ====="
python src/mosaic/export_community_usage.py "$PKL" "$ROOT/usage" "${SAMPLES[@]}"

echo "===== Stage G supplement (R) ====="
R_ENV
Rscript src/mosaic/mosaic_supplement.R "$ROOT/ssgsea" "$ROOT/usage" "${SAMPLES[@]}"

echo "PRODUCTION_OK: figure=$(ls -la results/fig_mosaic_nmf_supplement.pdf 2>/dev/null | awk '{print $5}') bytes"

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"mosaic_supplement (#40) $SLURM_JOB_ID finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
