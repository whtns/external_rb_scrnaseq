#!/bin/bash
#SBATCH --job-name=mosaic_smoke
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=03:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_smoke_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_smoke_%j.err

# Issue #40 SMOKE: run the WHOLE mosaicMPI chain (A->G) on TWO paper-retained
# samples with a tiny rank sweep (k=2..6, n_iter=20), on the debug partition, to
# validate correctness before the full 7-sample run:
#   - counts slot is raw/integer/symbol (export_counts.R asserts)
#   - a dense-.X mosaicMPI Dataset h5ad runs model-odg/factorize (no sparse error)
#   - check-h5ad passes
#   - integrate produces a network pkl and get_community_usage yields a barcode-
#     indexed matrix that joins 100% onto the Seurat cells (mosaic_supplement.R asserts)
#   - the supplement figure/table/methods are produced
# Two samples are used (not one) so `integrate` + community usage are exercised.

set -euo pipefail
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

SAMPLES=(SRX11133594 SRX11133593)
ROOT=output/mosaicmpi_smoke   # separate from the production output/mosaicmpi_rb
rm -rf "$ROOT"                # pure scratch -- start clean (avoids stale cNMF state)
mkdir -p logs "$ROOT/inputs"
GMT=data/mosaic_reference_programs.gmt
KFIRST=2; KLAST=4       # smoke: k=2,3,4 (production sweeps wider)
NITER=10                # smoke: few replicates (production ~100)
SEED=42
NGENES=2000             # top-N overdispersed genes -- default minscore kept 10.8k
                        # (way too many); ~2000 HVGs is standard cNMF and ~5x faster
NW=$SLURM_CPUS_PER_TASK  # parallel factorize workers

# r/4.4.1 lives under the gcc/13.3.0 Lmod hierarchy -- a bare `module purge` hides it
# ("exist but cannot be loaded as requested"), so reload gcc before r. Conda still gets
# a clean purge first (AGENTS.md: module-provided and conda python otherwise shadow).
R_ENV() {
  module purge
  module load gcc/13.3.0
  module load r/4.4.1
  export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}
MOSAIC_ENV() {
  module purge
  source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
  conda activate mosaic
  # One BLAS thread PER worker: we fan factorize out across $NW parallel processes,
  # so each cNMF process must be single-threaded or 16 procs x 16 BLAS threads
  # oversubscribes 16 cores (~16 min/task instead of ~1). This is the key to the
  # parallel speedup.
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
  mosaicmpi check-h5ad   -i "$ROOT/inputs/${s}.raw.h5ad" -o "$ROOT/inputs/${s}.h5ad"
  mosaicmpi model-odg    -n "$s" -o "$ROOT" -i "$ROOT/inputs/${s}.h5ad"
  # set-parameters via API (CLI's default_topn is broken in 2.5.0); k_first k_last:
  python src/mosaic/set_parameters.py "$ROOT" "$s" "$NGENES" $KFIRST $KLAST "$NITER" "$SEED"
  # Parallel factorize: launch NW workers, each running its slice of the k x iter
  # tasks. Serial factorize was ~2.5 min/task -> hours; this fans out across cpus.
  for i in $(seq 0 $((NW-1))); do
    mosaicmpi factorize -n "$s" -o "$ROOT" --worker_index "$i" --total_workers "$NW" &
  done
  wait
  mosaicmpi postprocess  -n "$s" -o "$ROOT" --cpus 1 --local_density_threshold 2.0
done

echo "===== Stage D create-config + integrate ====="
CFG_ARGS=(); for s in "${SAMPLES[@]}"; do CFG_ARGS+=(-i "$ROOT/$s/${s}.h5ad"); done
mosaicmpi create-config "${CFG_ARGS[@]}" -o "$ROOT/config.toml"
mosaicmpi integrate -o "$ROOT/integrate_all" -c "$ROOT/config.toml" --cpus "$SLURM_CPUS_PER_TASK"
PKL=$(ls "$ROOT/integrate_all"/network_integration.pkl* 2>/dev/null | head -1)
echo "integration pkl: $PKL"

echo "===== Stage E ssgsea ====="
mosaicmpi ssgsea -n "$PKL" -g "$GMT" -o "$ROOT/ssgsea" --cpus "$SLURM_CPUS_PER_TASK"

echo "===== Stage F community usage ====="
python src/mosaic/export_community_usage.py "$PKL" "$ROOT/usage" "${SAMPLES[@]}"

echo "===== Stage G supplement (R) ====="
R_ENV
Rscript src/mosaic/mosaic_supplement.R "$ROOT/ssgsea" "$ROOT/usage" "${SAMPLES[@]}"

echo "SMOKE_OK: figure=$(ls -la results/fig_mosaic_nmf_supplement.pdf 2>/dev/null | awk '{print $5}') bytes"
