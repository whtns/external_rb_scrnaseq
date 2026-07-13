#!/usr/bin/env Rscript
# Why does CellCycleScoring fail for these samples' filtered_seu?
# Report assays, default assay, cc-gene presence per assay, and the actual
# error from a direct CellCycleScoring attempt.
suppressPackageStartupMessages({ library(Seurat) })

files <- c(
  "output/seurat/SRX10831286_filtered_seu.rds",
  "output/seurat/SRX11133589_filtered_seu.rds"
)
sg <- cc.genes$s.genes; gg <- cc.genes$g2m.genes

for (f in files) {
  cat("\n================================================================\n", f, "\n")
  if (!file.exists(f)) { cat("  MISSING\n"); next }
  seu <- readRDS(f)
  cat("  n_cells:", ncol(seu), " assays:", paste(Assays(seu), collapse=","),
      " default:", DefaultAssay(seu), "\n")
  for (a in Assays(seu)) {
    rn <- rownames(seu[[a]])
    cat(sprintf("  assay %-6s: n_features=%d  s.genes present=%d/%d  g2m.genes present=%d/%d\n",
                a, length(rn), sum(sg %in% rn), length(sg), sum(gg %in% rn), length(gg)))
  }
  cat("  --- CellCycleScoring attempt (default assay) ---\n")
  res <- tryCatch({
    s2 <- Seurat::CellCycleScoring(seu, s.features = sg, g2m.features = gg, set.ident = FALSE)
    cat("  SUCCESS. S.Score range:", paste(round(range(s2$S.Score),3), collapse=".."), "\n")
    "ok"
  }, error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); "err" })

  if (identical(res, "err") && "gene" %in% Assays(seu)) {
    cat("  --- retry on 'gene' assay ---\n")
    tryCatch({
      DefaultAssay(seu) <- "gene"
      s3 <- Seurat::CellCycleScoring(seu, s.features = sg, g2m.features = gg, set.ident = FALSE)
      cat("  gene-assay SUCCESS. S.Score range:", paste(round(range(s3$S.Score),3), collapse=".."), "\n")
    }, error = function(e) cat("  gene-assay ERROR:", conditionMessage(e), "\n"))
  }
  rm(seu); gc(FALSE)
}
cat("\nDONE\n")
