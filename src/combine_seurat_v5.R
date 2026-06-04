#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)

one_q_sixteen_q_seus <- c(
	"SRX11133594" = "output/seurat/SRX11133594_filtered_seu.rds",
	"SRX11133593" = "output/seurat/SRX11133593_filtered_seu.rds",
	"SRX11133592" = "output/seurat/SRX11133592_filtered_seu.rds"
) %>% 
	map(readRDS)

scna_sets <- 
	list(
		"diploid_16q-" = c("", "16q-"),
		"16q-_16q-,1q+" = c("16q-", "16q- 1q+"),
		"diploid_16q-,1q+" = c("", "16q- 1q+")
	)

test1 <- map(one_q_sixteen_q_seus, subset, subset = scna %in% scna_sets[[3]])

test2 <- seuratTools::integration_workflow(test1)

# seu_list <- readRDS("output/seurat/g1_seu_list.rds")
#
# test0 <- purrr::reduce(seu_list, merge)

seu <- readRDS("output/seurat/g1_seu.rds")

DefaultAssay(seu) <- "gene"

seu$sample_id <- dplyr::coalesce(seu$sample_id, seu$sample_id.x)

seu <- seu[,seu$sample_id %in% c("SRX11133594", "SRX11133593", "SRX11133592")]

options(future.globals.maxSize = 3e+09)

seu[["RNA"]] <- seu[["gene"]]

DefaultAssay(seu) <- "RNA"

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample_id)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
seu <- FindClusters(seu, resolution = 2, cluster.name = "unintegrated_clusters")

seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

seu <- IntegrateLayers(
  seu, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:30)
seu <- FindClusters(seu, resolution = 2, cluster.name = "cca_clusters")

seu <- RunUMAP(seu, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca", min.dist = 0.05)

p1 <- DimPlot(
  seu,
  reduction = "umap.cca",
  group.by = c("sample_id", "clusters"),
  combine = FALSE, label.size = 2
)
