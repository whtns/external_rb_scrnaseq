#!/bin/bash
#SBATCH --job-name=mosaic_tail
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --time=00:40:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_tail_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/mosaic_tail_%j.err
set -euo pipefail
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
ROOT=output/mosaicmpi_smoke
GMT=data/mosaic_reference_programs.gmt
SAMPLES=(SRX11133594 SRX11133593)
PKL="$ROOT/integrate_all/network_integration.pkl.gz"

echo "===== Stage E ssgsea ====="
module purge; source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh; conda activate mosaic
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mosaicmpi ssgsea -n "$PKL" -g "$GMT" -o "$ROOT/ssgsea" --cpus "$SLURM_CPUS_PER_TASK"
echo "===== Stage F community usage ====="
python src/mosaic/export_community_usage.py "$PKL" "$ROOT/usage" "${SAMPLES[@]}"
echo "===== Stage G supplement ====="
module purge; module load gcc/13.3.0; module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
Rscript src/mosaic/mosaic_supplement.R "$ROOT/ssgsea" "$ROOT/usage" "${SAMPLES[@]}"
echo "TAIL_OK: figure=$(ls -la results/fig_mosaic_nmf_supplement.pdf 2>/dev/null | awk '{print $5}') bytes"
