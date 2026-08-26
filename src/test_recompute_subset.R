# Validate the FULL subset_seu_by_expression(recompute_pca=TRUE) path on one
# sample before the heavy split rebuild. Writes to an isolated test dir so it
# does not clobber the production low object; the seu_cells / hash DB rows are
# keyed by (temp) filepath so they don't collide with production either.
suppressPackageStartupMessages({
  source("packages.R"); library(Seurat); library(dplyr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

testdir <- "output/seurat/_recompute_test"
dir.create(testdir, showWarnings = FALSE, recursive = TRUE)
lab_src <- "output/seurat/SRX10264526_seu_hypoxia_labeled.rds"
lab     <- file.path(testdir, "SRX10264526_seu_hypoxia_labeled.rds")
stopifnot(file.exists(lab_src))
file.copy(lab_src, lab, overwrite = TRUE)

out <- subset_seu_by_expression(
  lab, run_hypoxia_clustering = TRUE,
  hypoxia_expr = "hypoxia_partition == 'low'", slug = "hypoxia_low",
  assay = "gene", recompute_pca = TRUE)
cat("\nwrote:", out, "\n")

s <- readRDS(out)
cat("cells:", ncol(s), "| reductions:", paste(names(s@reductions), collapse = ","), "\n")
cat("gene_snn_res.0.2 present (dictionary col):",
    "gene_snn_res.0.2" %in% colnames(s@meta.data), "\n")
cat("\n=== SCT clusters (fresh-PCA, read by collages) ===\n")
for (r in c(0.2, 0.4, 0.6, 0.8, 1)) {
  col <- paste0("SCT_snn_res.", r)
  if (col %in% colnames(s@meta.data)) {
    tab <- sort(table(as.character(s@meta.data[[col]])), decreasing = TRUE)
    cat(sprintf("  %s: %d clusters | min %d | median %.0f | n<20 = %d\n",
                col, length(tab), min(tab), median(tab), sum(tab < 20)))
  } else cat("  ", col, " MISSING\n")
}
cat("\nTEST DONE\n")
