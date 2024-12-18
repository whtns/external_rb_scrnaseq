#!/usr/bin/env Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
	library('tidyverse')
	library('fs')
	library(clustifyr)
	library(Seurat)
	library(seuratTools)
})


celltype_ref <- readRDS(celltype_ref)

counts <- Seurat::Read10X(matrix_dir)

seu <- Seurat::CreateSeuratObject(counts, assay = "gene")


DefaultAssay(seu) <- "gene"

seu <- 
	seu %>% 
	seurat_preprocess() %>% 
	seurat_reduce_dimensions()

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)

seu <- seurat_cluster(seu = seu, resolution = c(0.2, 0.4, 0.6),
											reduction = "pca")

seu <- SCTransform(seu, assay = "gene", verbose = FALSE)

seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)

DefaultAssay(seu) <- "SCT"

seu <- seurat_cluster(seu = seu, resolution = c(0.2, 0.4, 0.6),
                      reduction = "pca")

seu <- find_all_markers(seu, seurat_assay = "SCT")

seu <- add_percent_mito(seu, seurat_assay = "SCT")

# seu <- seuratTools::clustering_workflow(seu)

plot_celltype_predictions <- function(seu_path, celltype_ref, bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))) {
  # browser()

  seu_mat <- GetAssayData(seu, slot = "data")

  res <- clustify(
    input = seu_mat,
    metadata = seu$SCT_snn_res.0.2,
    ref_mat = celltype_ref,
    query_genes = VariableFeatures(seu)
  )

  cor_to_call(res)

  res2 <- cor_to_call(
    cor_mat = res,                  # matrix correlation coefficients
    cluster_col = "SCT_snn_res.0.2" # name of column in meta.data containing cell clusters
  )

  # Insert into original metadata as "type" column
  seu@meta.data <- call_to_metadata(
    res = res2,                     # data.frame of called cell type for each cluster
    metadata = seu@meta.data,           # original meta.data table containing cell clusters
    cluster_col = "SCT_snn_res.0.2" # name of column in meta.data containing cell clusters
  )

  return(seu)
}

seu <- plot_celltype_predictions(seu, celltype_ref)

# drop poor qc cells
# read count
seu <- seu[,seu$nCount_gene > 1000]
# genes detected
seu <- seu[,seu$nFeature_gene > 1000]
# %mito 
seu <- seu[,seu$percent.mt < 5]

saveRDS(seu, seu_path)


