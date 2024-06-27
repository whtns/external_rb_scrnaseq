#!/usr/bin/env Rscript

source("packages.R")
source("functions.R")
library(targets)

tar_load("filtered_seus")

pull_common_markers <- function(filtered_seus){
  names(filtered_seus) <- str_extract(filtered_seus, "SRR[0-9]*")

  # load filtered_seus ------------------------------
  # find cluster markers of every seu, compare with zinovyev markers
  # find common markers for clusters
  pull_cluster_markers <- function(seu_path){
    seu <- readRDS(seu_path)

    table_cluster_markers(seu)

  }

  my_cluster_markers <- map(filtered_seus, pull_cluster_markers)

  common_genes <-
    map(my_cluster_markers, "SCT_snn_res.0.4") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(abs(Average.Log.Fold.Change) > 0.5) %>%
    dplyr::arrange(Gene.Name) %>%
    dplyr::group_by(Gene.Name) %>%
    dplyr::filter(dplyr::n() > 3) %>%
    dplyr::filter(!str_detect(Gene.Name, pattern = "^RP.*")) %>%
    dplyr::filter(!str_detect(Gene.Name, pattern = "^MT.*")) %>%
    dplyr::distinct(Gene.Name) %>%
    dplyr::ungroup() %>%
    # dplyr::slice_sample(n =50) %>%
    # dplyr::pull(Gene.Name) %>%
    # unique() %>%
    # sample(50) %>%
    identity()

}

common_genes <- pull_common_markers(filtered_seus)

heatmap_marker_genes <- function(seu_path, common_genes){
  # browser()
  sample_id <- str_extract(seu_path, "SRR[0-9]*")
  seu <- readRDS(seu_path)

  test0 <- seu@misc$markers$SCT_snn_res.0.4$presto %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_head(n=10) %>%
    dplyr::select(Gene.Name, Cluster) %>%
    dplyr::inner_join(common_genes, by = "Gene.Name") %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Cluster) %>%
    identity()

  mymarkers <-
    test0 %>%
    dplyr::pull(Gene.Name)

  seu$scna <- na_if(seu$scna, "")
  seu$scna <- replace_na(seu$scna, "none")

  seu_heatmap <- ggplotify::as.ggplot(seu_complex_heatmap(seu, features = mymarkers, group.by = c("SCT_snn_res.0.4", "zphase", "scna"), col_arrangement = c("SCT_snn_res.0.4", "zphase", "scna"), cluster_rows = FALSE)) +
    labs(title = sample_id)

  # seu_heatmap <- Seurat::DoHeatmap(seu, features = mymarkers, group.by = "SCT_snn_res.0.4")

  # ggsave(glue("results/{sample_id}_heatmap.pdf"), height = 8, width = 8)

  return(list("markers" = test0, "heatmap" = seu_heatmap))

}

# test1 <- heatmap_marker_genes(filtered_seus[[1]], common_genes)

seurat_heatmaps <- map(filtered_seus, heatmap_marker_genes, common_genes)

test0 <- map(seurat_heatmaps, "markers") %>%
  map(tidyr::pivot_wider, names_from = "Cluster", values_from = "Gene.Name")

pdf("results/seurat_heatmaps.pdf")
seurat_heatmaps %>%
  map("heatmap")
dev.off()

# zinovyev_genes <- read_tsv("data/zinovyev_cc_genes.tsv")
#
# plotted_genes <- c(common_genes, zinovyev_genes$symbol)
