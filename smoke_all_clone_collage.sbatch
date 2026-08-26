#!/bin/bash
#SBATCH --job-name=smoke_allclone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:55:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_allclone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_allclone_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_all_clone_collage.R
