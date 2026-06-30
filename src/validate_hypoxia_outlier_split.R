#!/usr/bin/env Rscript
# Validate the outlier-mean-score hypoxia split rule on the user's reference
# samples. Replicates the single-round multi-resolution clustering done by
# split_hypoxia_by_clusters() and prints the per-cluster decision table from the
# updated identify_hypoxia_clusters(), WITHOUT the expensive subset/recluster
# step. Confirms the expected clusters flag as high-side outliers at res 0.2.
suppressPackageStartupMessages(source("packages.R"))
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)

# sample -> expected is_outlier clusters at resolution 0.2 (user-supplied)
expected <- list(
  SRX10264517 = c("2", "4"),
  SRX10264518 = c("4", "5"),
  SRX10264519 = c("2"),
  SRX10264520 = c("2")
)

resolutions <- c(0.2, 0.4, 0.6)
split_assay <- "gene"
out_dir <- "results/hypoxia_cluster_split"
fs::dir_create(out_dir)

for (sid in names(expected)) {
  seu_path <- file.path("output/seurat", paste0(sid, "_seu_hypoxia.rds"))
  cat("\n=================== ", sid, " ===================\n")
  if (!file.exists(seu_path)) { cat("MISSING input:", seu_path, "\n"); next }

  seu <- readRDS(seu_path)
  da <- if (split_assay %in% names(seu@assays)) split_assay else Seurat::DefaultAssay(seu)
  Seurat::DefaultAssay(seu) <- da
  if (!"data" %in% SeuratObject::Layers(seu[[da]]))
    seu <- Seurat::NormalizeData(seu, assay = da, verbose = FALSE)

  if (!"pca" %in% names(seu@reductions)) { cat("no pca reduction; skip\n"); next }

  n_cells <- ncol(seu)
  graph_names <- paste0(da, c("_nn", "_snn"))
  seu <- Seurat::FindNeighbors(seu, dims = 1:30, reduction = "pca",
                               graph.name = graph_names,
                               k.param = min(20L, n_cells - 1L), verbose = FALSE)

  det_all <- list()
  for (res in resolutions) {
    grp <- glue::glue("hypsplit_r1_res.{res}")
    seu <- Seurat::FindClusters(seu, graph.name = paste0(da, "_snn"),
                                resolution = res, verbose = FALSE)
    seu@meta.data[[grp]] <- seu$seurat_clusters
    det <- identify_hypoxia_clusters(seu, grp, return_detail = TRUE)
    det <- cbind(sample_id = sid, round = 1, resolution = res, pool_n = n_cells,
                 det, stringsAsFactors = FALSE)
    det_all[[as.character(res)]] <- det
  }
  det_all <- dplyr::bind_rows(det_all)
  readr::write_csv(det_all, file.path(out_dir, paste0(sid, "_validate_split_log.csv")))

  # focus on res 0.2
  d02 <- det_all[det_all$resolution == 0.2, ]
  cat("--- resolution 0.2 ---\n")
  print(d02[, c("cluster", "n_cells", "top_markers", "n_hyp_markers",
                "mean_score", "outlier_fence", "is_outlier")], row.names = FALSE)
  got <- sort(as.character(d02$cluster[d02$is_outlier]))
  exp <- sort(expected[[sid]])
  ok <- identical(got, exp)
  cat(sprintf("\nEXPECTED outliers @res0.2: {%s}\n", paste(exp, collapse = ", ")))
  cat(sprintf("GOT      outliers @res0.2: {%s}   => %s\n",
              paste(got, collapse = ", "), if (ok) "MATCH" else "MISMATCH"))
}
cat("\n==== DONE ====\n")
