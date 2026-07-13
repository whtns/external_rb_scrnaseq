#!/usr/bin/env Rscript
# Render, per Collin sample, a cluster UMAP DimPlot beside a table of the top-10
# marker genes per cluster. Cluster labels and the marker table both come from
# gene_snn_res.0.2 (the resolution the cluster_markers DB was computed at), so the
# DimPlot clusters and the table rows correspond 1:1.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
  library(dplyr); library(readr)
})

samples <- c("SRX10031191", "SRX10031192", "SRX10031194")
outdir  <- "results/cluster_markers"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (s in samples) {
  seu_path <- sprintf("output/seurat/%s_seu.rds", s)
  csv_path <- sprintf("%s/%s_cluster_markers.csv", outdir, s)
  message("== ", s, " ==")
  seu <- readRDS(seu_path)
  mk  <- readr::read_csv(csv_path, show_col_types = FALSE)

  # Choose the resolution column that matches the DB cluster count (prefer 0.2).
  db_nclus <- length(unique(mk$cluster))
  res_cols <- grep("_snn_res", colnames(seu@meta.data), value = TRUE)
  clus_col <- if ("gene_snn_res.0.2" %in% res_cols) {
    "gene_snn_res.0.2"
  } else {
    counts <- vapply(res_cols, function(cc) length(unique(seu@meta.data[[cc]])), integer(1))
    hit <- res_cols[which(counts == db_nclus)]
    if (length(hit)) hit[[1]] else res_cols[[1]]
  }
  message("  cluster column: ", clus_col,
          " (seu clusters=", length(unique(seu@meta.data[[clus_col]])),
          ", db clusters=", db_nclus, ")")

  lvls <- as.character(sort(unique(as.integer(as.character(seu@meta.data[[clus_col]])))))
  seu$clusters <- factor(as.character(seu@meta.data[[clus_col]]), levels = lvls)

  red <- if ("umap" %in% names(seu@reductions)) "umap" else names(seu@reductions)[[1]]
  dplot <- DimPlot(seu, group.by = "clusters", reduction = red, label = TRUE, repel = TRUE) +
    labs(title = sprintf("%s — clusters (%s, %s)", s, clus_col, red)) +
    theme(legend.position = "right")

  # Marker table rendered as a scalable ggplot text grid (cluster x rank).
  tab <- mk %>% dplyr::filter(marker_rank <= 10)
  gtab <- ggplot(tab, aes(x = factor(marker_rank),
                          y = factor(cluster, levels = rev(lvls)))) +
    geom_text(aes(label = gene_name), size = 2.6) +
    scale_x_discrete(position = "top") +
    labs(x = "marker rank", y = "cluster",
         title = sprintf("%s — top-10 cluster markers (%s)", s, clus_col)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold"))

  combined <- (dplot | gtab) + patchwork::plot_layout(widths = c(1, 1.2))
  out <- sprintf("%s/%s_clusters_and_markers.pdf", outdir, s)
  ggsave(out, combined, width = 16, height = 8)
  cat("wrote ", out, "\n", sep = "")
}
cat("done\n")
