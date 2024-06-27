#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(chevreul)

regressed_seus <- dir_ls("output/seurat", glob = "*regressed_seu.rds", recurse = TRUE) %>%
  set_names(str_extract(., "SRR[0-9]*")) %>%
  map(readRDS)

saveRDS(seurat_integration_pipeline(regressed_seus), "output/seurat/seurat_merged_seu.rds")


combined_seu <- readRDS("output/seurat/seurat_merged_seu.rds")

regressed_combined_seu <- combined_seu

regressed_combined_seu@misc$markers$gene_snn_res.0.2 <- NULL

regressed_combined_seu <- CellCycleScoring(regressed_combined_seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

regressed_combined_seu <- ScaleData(regressed_combined_seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(regressed_combined_seu))
regressed_combined_seu <- seurat_reduce_dimensions(regressed_combined_seu)

regressed_combined_seu <- seurat_cluster(seu = regressed_combined_seu, resolution = 0.2,
                      reduction = "pca")

regressed_combined_seu <- find_all_markers(regressed_combined_seu, seurat_assay = "integrated")

plot_scna_per_sample <- function(sample_id, seu){
  sample_seu <- seu[,seu$batch == sample_id]
  DimPlot(sample_seu, group.by = c("scna"), combine = FALSE)
}

sample_ids <-
  unique(combined_seu$batch) %>%
  set_names(.)

scna_plots <- map(sample_ids, plot_scna_per_sample, regressed_combined_seu)

pdf("results/seu_combined_scna_plots.pdf")
scna_plots
dev.off()

browseURL("results/seu_combined_scna_plots.pdf")

plot_markers(regressed_combined_seu, "integrated_snn_res.0.2", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10, unique_markers = TRUE)

ggsave("results/integrated_regressed_markers.pdf", height = 12, width = 8)

browseURL("results/integrated_regressed_markers.pdf")

DimPlot(regressed_combined_seu, group.by = "integrated_snn_res.0.2", split.by = "Phase")
ggsave("results/integrated_regressed_louvain_split_phase.pdf", height = 8, width = 12)
browseURL("results/integrated_regressed_louvain_split_phase.pdf")


