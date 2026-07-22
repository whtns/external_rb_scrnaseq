#!/bin/bash
#SBATCH --job-name=derive_dipchroms
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/derive_dipchroms_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/derive_dipchroms_%j.err

# Derive a per-sample diploid_chroms set from each sample's own call burden.
# Read-only over output/numbat_sridhar; writes
# results/diploid_chroms_per_sample.tsv, which rerun_numbat_diploidpin.sbatch
# consumes. Run this BEFORE that array job.
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/derive_diploid_chroms.R
