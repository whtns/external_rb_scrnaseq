#!/bin/bash
#SBATCH --job-name=backfill_cells_db
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=60:00
#SBATCH --output=logs/backfill_cells_db_%j.out
#SBATCH --error=logs/backfill_cells_db_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript backfill_cell_barcodes_db.R
