#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(chevreul)
library(patchwork)
library(stringr)
library(glue)



seu_paths <- dir_ls("output/seurat/", glob = "*filtered_seu.rds", recurse = TRUE)

plot_effects_of_cell_cycle_regression <- function(seu_path) {

  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  regressed_seu_path <- str_replace(seu_path, "_filtered_seu.rds", "_regressed_seu.rds")
  filtered_seu <- readRDS(seu_path)

  # filtered_seu0 <- filtered_seu
  # filtered_seu0@misc$markers$gene_snn_res.0.2 <- NULL
  #
  # filtered_seu0 <- ScaleData(filtered_seu0, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(filtered_seu0))
  # filtered_seu0 <- seurat_reduce_dimensions(filtered_seu0)
  #
  # filtered_seu0 <- seurat_cluster(seu = filtered_seu0, resolution = 0.2,
  #                       reduction = "pca")
  #
  # filtered_seu0 <- find_all_markers(filtered_seu0, seurat_assay = "gene")
  #
  # saveRDS(filtered_seu0, regressed_seu_path)

  filtered_seu0 <- readRDS(regressed_seu_path)

  fs::dir_create("results/effect_of_regression")

  # plot_markers(filtered_seu, "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
  #   labs(title = "original") +
  # plot_markers(filtered_seu0, "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_markers.pdf"), height = 10, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_markers.pdf"))

  filtered_seu$scna[filtered_seu$scna == ""] <- "none"
  filtered_seu0$scna[filtered_seu0$scna == ""] <- "none"

  plot_markers(filtered_seu, "scna", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10, unique_markers = TRUE) +
    labs(title = sample_id)

  ggsave(glue("results/effect_of_regression/{sample_id}_markers_by_scna.pdf"), height = 8, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_markers_by_scna.pdf"))
  #
  # #   original_mt_plot <-
  # #     FeaturePlot(filtered_seu, features = "percent.mt") +
  # #     labs(title = "original")
  # #
  # #   regressed_mt_plot <- FeaturePlot(filtered_seu0, features = "percent.mt") +
  # #   labs(title = "regressed")
  # #
  # #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  # #     plot_annotation(title = sample_id)
  # #
  # # ggsave(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"), heigh = 6, width = 10)
  # # browseURL(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"))
  #
  # filtered_seu <- AddModuleScore(filtered_seu, subtype_markers)
  # filtered_seu0 <- AddModuleScore(filtered_seu0, subtype_markers)
  #
  #   original_mt_plot <-
  #     FeaturePlot(filtered_seu, features = "Cluster1") +
  #     labs(title = "original")
  #
  #   regressed_mt_plot <- FeaturePlot(filtered_seu0, features = "Cluster1") +
  #   labs(title = "regressed")
  #
  #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #     plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"))
  #
  # original_mt_plot <-
  #   FeaturePlot(filtered_seu, features = "Cluster2") +
  #   labs(title = "original")
  #
  # regressed_mt_plot <- FeaturePlot(filtered_seu0, features = "Cluster2") +
  #   labs(title = "regressed")
  #
  # wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"))

  # DimPlot(filtered_seu, group.by = "gene_snn_res.0.2") +
  #   labs(title = "original") +
  #   DimPlot(filtered_seu0, group.by = "gene_snn_res.0.2") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_louvain.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_louvain.pdf"))
  #
  # DimPlot(filtered_seu, group.by = "scna") +
  #   labs(title = "original") +
  #   DimPlot(filtered_seu0, group.by = "scna") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_scna.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_scna.pdf"))
  #
  # DimPlot(filtered_seu, group.by = "Phase") +
  #   labs(title = "original") +
  #   DimPlot(filtered_seu0, group.by = "Phase") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_phase.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_phase.pdf"))

  return(regressed_seu_path)

}

map(seu_paths, plot_effects_of_cell_cycle_regression)

# plot_effects_of_cell_cycle_regression(seu_paths[9])
