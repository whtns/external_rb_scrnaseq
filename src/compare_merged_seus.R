#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(Seurat)
library(chevreul)
library(sctransform)

seu_paths <- c(
  "output/seurat/SRR14800534_regressed_seu.rds",
  "output/seurat/SRR14800535_regressed_seu.rds",
  "output/seurat/SRR14800536_regressed_seu.rds"
  )

sct_prep <- function(seu_path){
  seu <- readRDS(seu_path)

  seu <- SCTransform(seu, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.7, verbose = FALSE)
}

seu_list <- map(seu_paths, sct_prep)

features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = features)

seu.anchors <- FindIntegrationAnchors(object.list = seu_list, normalization.method = "SCT",
                                         anchor.features = features)
seu.combined.sct <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT")

seu.combined.sct <- RunPCA(seu.combined.sct, verbose = FALSE)
seu.combined.sct <- RunUMAP(seu.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
seu.combined.sct <- FindNeighbors(seu.combined.sct, reduction = "pca", dims = 1:30)
seu.combined.sct <- FindClusters(seu.combined.sct, resolution = 0.3)





merged_seu <- readRDS("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds")
