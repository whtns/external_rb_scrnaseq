#!/usr/bin/env Rscript
# Which score columns do the FILTERED objects (the ones the two-clone filtered
# collages are built from) actually carry? Decides whether
# plot_scna_two_clone_res_collages has to recompute hypoxia / MT before it can
# annotate the heatmap columns with them.
suppressPackageStartupMessages({
  library(Seurat); library(stringr)
})

want <- c("hypoxia", "MT", "hypoxia_score", "percent.mt", "S.Score", "G2M.Score")

paths <- c(
  "output/seurat/SRX10264523_filtered_seu.rds",
  "output/seurat/SRX11133594_filtered_seu.rds",
  "output/seurat/SRX10264523_hypoxia_low_seu.rds"
)

for (p in paths) {
  if (!file.exists(p)) { cat(p, "| NO FILE\n"); next }
  seu <- readRDS(p)
  md  <- seu@meta.data
  present <- want[want %in% colnames(md)]
  missing <- setdiff(want, present)
  cat("\n==", basename(p), "| ncell =", ncol(seu), "\n")
  cat("   present:", paste(present, collapse = ", "), "\n")
  cat("   MISSING:", paste(missing, collapse = ", "), "\n")
  for (cl in present) {
    v <- md[[cl]]
    if (is.numeric(v))
      cat(sprintf("   %-14s range %.4f .. %.4f  (NA %d)\n", cl,
                  min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(is.na(v))))
  }
  rm(seu); gc()
}
cat("\nDIAG DONE\n")
