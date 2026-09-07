#!/bin/bash
#SBATCH --job-name=diag_dosage
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/diag_clone_dosage_depth_%j.log

# #43 follow-up part 2: does genome dosage explain the depth-clone association?
# Read-only -- reads output/numbat_sridhar/*_numbat.rds and the cell_qc_values
# table, writes only results/. Touches no numbat outputs, so the SRX11133592/93/94
# hold-out is unaffected. Does NOT touch the _targets_r431 store, so it is safe
# to run alongside a pipeline job.
#
# 64G because an earlier attempt in a 1-cpu interactive allocation was OOM-killed
# (exit 137) three samples in; part 1 of this analysis peaked at 4.5 GB under
# the same 64G ask.

set -euo pipefail

module purge
module load gcc/13.3.0   # r/4.4.1 sits under gcc in the Lmod hierarchy;
module load r/4.4.1      # a bare purge hides it (job 11784414 failed on this)

export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

echo "=== stage 1: smoke on two samples ==="
Rscript src/diag_clone_dosage_depth.R SRX10264519 SRX22868104
echo "=== smoke OK ==="

echo "=== stage 2: full cohort ==="
Rscript src/diag_clone_dosage_depth.R
echo "=== done ==="
