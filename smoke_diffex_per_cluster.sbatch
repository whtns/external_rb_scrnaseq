#!/bin/bash
#SBATCH --job-name=smoke_diffex
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_diffex_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/smoke_diffex_%j.err

# Read-only smoke of the diffex-per-cluster function on one sample (load_all, no
# tar_make) -- safe to run alongside another tar_make job.
module load r/4.4.1
module load poppler/23.04.0
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/smoke_diffex_per_cluster.R
