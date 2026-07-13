#!/bin/bash
#SBATCH --job-name=test_hyp_split
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/test_hyp_split_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/test_hyp_split_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/test_hypoxia_cluster_split.R "${1:-output/seurat/SRX10031194_seu_hypoxia.rds}"
