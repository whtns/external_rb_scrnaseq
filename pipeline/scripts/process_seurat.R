#!/usr/bin/env Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
	library('tidyverse')
	library('fs')
	library(Seurat)
	library(presto)
})

# ---------------------------------------------------------------------------
# Inlined seuratTools helpers (avoids installing the full package)
# ---------------------------------------------------------------------------

seurat_preprocess <- function(seu) {
	seu <- NormalizeData(seu, verbose = FALSE)
	seu <- FindVariableFeatures(seu, selection.method = "vst", verbose = FALSE)
	seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
	seu
}

seurat_reduce_dimensions <- function(seu) {
	npcs <- if (ncol(seu) < 50) ncol(seu) - 1L else 50L
	seu <- RunPCA(seu, assay = "gene", features = VariableFeatures(seu), npcs = npcs, verbose = FALSE)
	if ((ncol(seu) - 1) > 3 * 30) {
		seu <- RunTSNE(seu, assay = "gene", reduction = "pca", dims = 1:30, check_duplicates = FALSE)
		seu <- RunUMAP(seu, assay = "gene", reduction = "pca", dims = 1:30, verbose = FALSE)
	}
	seu
}

seurat_cluster <- function(seu, resolution = c(0.2, 0.4, 0.6), reduction = "pca") {
	seu <- FindNeighbors(seu, dims = 1:30, reduction = reduction, verbose = FALSE)
	for (r in resolution) {
		seu <- FindClusters(seu, resolution = r, algorithm = 1, verbose = FALSE)
	}
	seu
}

stash_marker_features <- function(metavar, seu, seurat_assay = "gene",
                                   top_n = 200, p_val_cutoff = 0.5) {
	message(paste0("stashing presto markers for ", metavar))
	Idents(seu) <- seu@meta.data[[metavar]]
	DefaultAssay(seu) <- seurat_assay
	list(
		presto = presto::wilcoxauc(seu, metavar, seurat_assay = seurat_assay) %>%
			dplyr::group_by(group) %>%
			dplyr::filter(padj < p_val_cutoff) %>%
			dplyr::top_n(n = top_n, wt = logFC) %>%
			dplyr::arrange(group, dplyr::desc(logFC)) %>%
			dplyr::select(Gene.Name = feature, Average.Log.Fold.Change = logFC,
			              Adjusted.pvalue = padj, avgExpr, Cluster = group)
	)
}

find_all_markers <- function(seu, seurat_assay = "gene") {
	if (inherits(seu[[seurat_assay]], "Assay5")) {
		seu[[seurat_assay]] <- JoinLayers(seu[[seurat_assay]])
	}
	metavar <- colnames(seu@meta.data)[grepl(paste0(seurat_assay, "_snn_res\\."), colnames(seu@meta.data))]
	new_markers <- purrr::map(metavar, stash_marker_features, seu = seu, seurat_assay = seurat_assay)
	names(new_markers) <- metavar
	seu@misc$markers <- c(seu@misc$markers[!names(seu@misc$markers) %in% metavar], new_markers)
	seu
}

add_percent_mito <- function(seu, seurat_assay = "gene") {
	seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-", assay = seurat_assay)
	seu
}

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

celltype_ref <- readRDS(celltype_ref)

counts <- Seurat::Read10X(matrix_dir)

seu <- Seurat::CreateSeuratObject(counts, assay = "gene")

DefaultAssay(seu) <- "gene"

seu <- add_percent_mito(seu, seurat_assay = "gene")
seu <- seu[, seu$nCount_gene > 1000]
seu <- seu[, seu$nFeature_gene > 1000]
seu <- seu[, seu$percent.mt < 10]

seu <-
	seu %>%
	seurat_preprocess() %>%
	seurat_reduce_dimensions()

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)

seu <- seurat_cluster(seu = seu, resolution = c(0.2, 0.4, 0.6),
											reduction = "pca")

options(future.globals.maxSize = Inf)
seu <- SCTransform(seu, assay = "gene", verbose = FALSE)

seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)

DefaultAssay(seu) <- "gene"

seu <- seurat_cluster(seu = seu, resolution = c(0.2, 0.4, 0.6),
                      reduction = "pca")

seu <- find_all_markers(seu, seurat_assay = "gene")

saveRDS(seu, seu_path)
