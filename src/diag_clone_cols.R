#!/usr/bin/env Rscript
# What clone-related metadata columns exist in a hypoxia_low seu, and their values?
suppressPackageStartupMessages({ library(Seurat) })
f <- "output/seurat/SRX11133594_hypoxia_low_seu.rds"
seu <- readRDS(f)
md <- seu@meta.data
cat(f, " n_cells:", ncol(seu), "\n")
cat("metadata columns:\n"); print(colnames(md))
for (col in c("clone_opt","clone","GT_opt","scna")) {
  if (col %in% colnames(md)) {
    cat("\n== ", col, " ==\n"); print(table(md[[col]], useNA="always"))
  } else cat("\n== ", col, ": ABSENT ==\n")
}
cat("\nDONE\n")
