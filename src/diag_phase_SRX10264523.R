#!/usr/bin/env Rscript
# Does SRX10264523_hypoxia_low have a real Seurat Phase, or all-G1 from
# zero cell-cycle scores? Report per group.by cluster: modal Phase and the
# Phase distribution, plus score ranges.
suppressPackageStartupMessages({ library(Seurat) })

f <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
seu <- readRDS(f)
md <- seu@meta.data
cat(f, "  n_cells:", ncol(seu), "\n")
cat("S.Score:", paste(round(range(md$S.Score, na.rm=TRUE),3), collapse=".."),
    " G2M.Score:", paste(round(range(md$G2M.Score, na.rm=TRUE),3), collapse=".."), "\n")
cat("Phase overall:\n"); print(table(md$Phase, useNA="always"))

grp <- grep("SCT_snn_res", colnames(md), value=TRUE)
cat("\nresolution cols:", paste(grp, collapse=", "), "\n")
for (g in c("SCT_snn_res.0.6", grp[1])) {
  if (!g %in% colnames(md)) next
  cat("\n=== modal Phase per", g, "cluster ===\n")
  tab <- table(md[[g]], md$Phase)
  modal <- apply(tab, 1, function(r) colnames(tab)[which.max(r)])
  for (cl in rownames(tab)) {
    cat(sprintf("  cluster %-4s modal=%-4s  dist: %s\n", cl, modal[cl],
                paste(sprintf("%s=%d", colnames(tab), tab[cl,]), collapse=" ")))
  }
  break
}
cat("\nDONE\n")
