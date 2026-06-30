#!/usr/bin/env Rscript
# Smoke test for split_hypoxia_by_clusters() on a single sample.
# Load the pipeline's full package set (provides seuratTools::find_all_markers,
# tidyverse, Seurat, etc.), then override numbatHelpers with the dev source so
# the new split_hypoxia_by_clusters() is picked up.
suppressPackageStartupMessages(source("packages.R"))
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
seu_path <- if (length(args) >= 1) args[[1]] else "output/seurat/SRX10031194_seu_hypoxia.rds"

cat("==== testing split on:", seu_path, "====\n")
out <- split_hypoxia_by_clusters(
  seu_path,
  resolutions = c(0.2, 0.4, 0.6),
  n_iter = 1,
  split_assay = "gene", low_assay = "gene", high_assay = "SCT",
  write_qc = TRUE
)

cat("\n==== returned paths ====\n")
print(out)

cat("\n==== labelled object summary ====\n")
labeled <- sub("_seu_hypoxia\\.rds$", "_seu_hypoxia_labeled.rds", seu_path)
seu <- readRDS(labeled)
cat("hypoxia_partition:\n"); print(table(seu$hypoxia_partition, useNA = "always"))
cat("hypoxia_round:\n"); print(table(seu$hypoxia_round, useNA = "always"))

cat("\n==== low / high cell counts ====\n")
for (nm in names(out)) {
  if (!is.na(out[[nm]]) && file.exists(out[[nm]])) {
    s <- readRDS(out[[nm]])
    cat(sprintf("%-5s %-55s n=%d\n", nm, basename(out[[nm]]), ncol(s)))
  } else {
    cat(sprintf("%-5s MISSING (%s)\n", nm, out[[nm]]))
  }
}

cat("\n==== top markers of high subset (sanity: expect hypoxia genes) ====\n")
if (!is.na(out[["high"]]) && file.exists(out[["high"]])) {
  sh <- readRDS(out[["high"]])
  DefaultAssay(sh) <- "gene"
  grp_cols <- grep("_snn_res\\.", colnames(sh@meta.data), value = TRUE)
  if (ncol(sh) >= 20 && length(grp_cols) > 0) {
    grp <- grp_cols[[1]]
    m <- tryCatch(presto::wilcoxauc(sh, grp, seurat_assay = "gene"), error = function(e) NULL)
    if (!is.null(m)) {
      top <- m |>
        dplyr::filter(padj < 0.5) |>
        dplyr::group_by(group) |>
        dplyr::arrange(desc(logFC), .by_group = TRUE) |>
        dplyr::slice_head(n = 8) |>
        dplyr::summarise(genes = paste(feature, collapse = ", "), .groups = "drop")
      print(as.data.frame(top))
    }
  } else {
    cat("high subset too small or no cluster columns for marker check\n")
  }
}
cat("\n==== DONE ====\n")
