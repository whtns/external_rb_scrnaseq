#!/bin/bash
#SBATCH --job-name=verify_hypoxia_counts
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/verify_hypoxia_counts_%j.log

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

MSG=$(Rscript -e '
suppressPackageStartupMessages(library(SeuratObject))
f <- function(p){ if(length(p)==0 || !file.exists(p)) return(NA_integer_); ncol(readRDS(p)) }
lo <- Sys.glob("output/seurat/*SRX11133593*hypoxia_low*.rds")[1]
hi <- Sys.glob("output/seurat/*SRX11133593*hypoxia_high*.rds")[1]
fl <- "output/seurat/SRX11133593_filtered_seu.rds"
nf <- f(fl); nl <- f(lo); nh <- f(hi)
cat(sprintf("SRX11133593 filtered=%d low=%d high=%d (low+high=%d, diff=%d)",
            nf, nl, nh, nl+nh, nf-(nl+nh)))
' 2>/dev/null)

echo "$MSG"
curl -s -X POST -H 'Content-type: application/json' \
  --data "{\"text\":\"hypoxia-chain rebuild done. $MSG\"}" "$(cat ~/.slack_webhook)"
