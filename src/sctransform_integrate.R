#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)

source("packages.R")
source("functions.R")

seu_paths <- c("output/seurat/SRR14800534_filtered_seu.rds",
               "output/seurat/SRR14800535_filtered_seu.rds",
               "output/seurat/SRR14800536_filtered_seu.rds"
) %>%
  set_names(str_extract(., "SRR[0-9]*"))

prep_seus <- function(seu_path){
  seu <-
    seu_path %>%
    readRDS() %>%
    SCTransform(assay = "gene", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.4, verbose = FALSE)

  return(seu)
}

# seus <- map(seu_paths, prep_seus)

seus <- map(seu_paths, readRDS) %>%
  imap(~{.x$batch <- .y; return(.x)})

features <- SelectIntegrationFeatures(object.list = seus, nfeatures = 3000)
seus <- PrepSCTIntegration(object.list = seus, anchor.features = features)

seu.anchors <- FindIntegrationAnchors(object.list = seus, normalization.method = "SCT",
                                         anchor.features = features)
seu.combined.sct <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT")


seu.combined.sct <- RunPCA(seu.combined.sct, verbose = FALSE)
seu.combined.sct <- RunUMAP(seu.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
seu.combined.sct <- FindNeighbors(seu.combined.sct, reduction = "pca", dims = 1:30)
seu.combined.sct <- FindClusters(seu.combined.sct, resolution = c(0.1, 0.2, 0.4))

seu.combined.sct <- find_all_markers(seu.combined.sct, "integrated_snn_res.0.2", seurat_assay = "integrated")

plot_markers(seu.combined.sct, metavar = "integrated_snn_res.0.2", num_markers = 10, return_plotly = FALSE)

seu.combined.sct <- find_all_markers(seu.combined.sct, "scna", seurat_assay = "integrated")
plot_markers(seu.combined.sct, "scna", num_markers = 10)

# plot_distribution_of_clones_across_clusters(seu.combined.sct, "asdf", var_x = "integrated_snn_res.0.2")

seu.combined.sct <- assign_zphase(seu.combined.sct)

saveRDS(seu.combined.sct, "output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds")

