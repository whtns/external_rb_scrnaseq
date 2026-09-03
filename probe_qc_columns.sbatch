#!/bin/bash
#SBATCH --job-name=probe_qccols
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/probe_qccols_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/probe_qccols_%j.err

# SRX10264520/21/22 ingest without error but write ZERO cell_qc_values rows,
# while SRX10031193 wrote 3046. extract_seu_metadata() populates that table from
# nCount_gene / nFeature_gene / percent.mt, so the likely cause is that these
# objects carry different metadata column names (e.g. only nCount_RNA).
#
# Read-only: reports the metadata columns, writes nothing.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
suppressPackageStartupMessages(library(Seurat))
for (s in c("SRX10264520", "SRX10031193")) {
  f <- file.path("output/seurat", paste0(s, "_seu.rds"))
  seu <- readRDS(f)
  md <- seu@meta.data
  cat("\n=== ", s, " (", nrow(md), " cells) ===\n", sep = "")
  cat("assays: ", paste(SeuratObject::Assays(seu), collapse = ", "), "\n", sep = "")
  hits <- grep("nCount|nFeature|percent", colnames(md), value = TRUE)
  cat("count/feature/percent columns: ", paste(hits, collapse = ", "), "\n", sep = "")
  for (c in hits) {
    v <- md[[c]]
    if (is.numeric(v)) cat(sprintf("   %-22s median=%.1f  n_NA=%d\n", c, median(v, na.rm=TRUE), sum(is.na(v))))
  }
  rm(seu, md); gc(verbose = FALSE)
}'
