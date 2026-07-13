#!/usr/bin/env Rscript
# Test auto_phase_level() on SRX10264523_hypoxia_low at SCT_snn_res.0.6.
suppressPackageStartupMessages({
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")
  library(Seurat)
})

f <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
seu <- readRDS(f)
g <- "SCT_snn_res.0.6"

lab <- auto_phase_level(seu, g)
cat("=== auto_phase_level (", g, ") ===\n")
md <- seu@meta.data
for (cl in names(lab)) {
  n <- sum(as.character(md[[g]]) == cl)
  ms <- mean(md$S.Score[as.character(md[[g]]) == cl], na.rm = TRUE)
  mg <- mean(md$G2M.Score[as.character(md[[g]]) == cl], na.rm = TRUE)
  cat(sprintf("  cluster %-3s -> %-7s  (n=%d, meanS=%.3f, meanG2M=%.3f)  =>  %s_%s\n",
              cl, lab[cl], n, ms, mg, lab[cl], cl))
}
cat("\nphase label tally:\n"); print(table(lab))
cat("\nDONE\n")
