#!/usr/bin/env Rscript 

source("packages.R")
source("functions.R")
library(targets)
library(SingleCellExperiment)
library(ComplexHeatmap)

tar_load("debranched_seus_1q")
tar_load("debranched_seus_2p")
tar_load("debranched_seus_6p")
tar_load("debranched_seus_16q")




read_seu_by_phase <- function(seu_path, phase_level = "g1"){
	# browser()
	seu <- readRDS(seu_path)
	
	DefaultAssay(seu) <- "gene"
	
	seu <- DietSeurat(
		seu,
		assays = "gene",
		dimreducs = c("pca"),
	)
	
	meta <- seu@meta.data
	
	if(is(seu$gene, "Assay")){
		seu_v5 <- CreateSeuratObject(counts = seu$gene@counts, 
																 data = seu$gene@data, assay = "gene", meta.data = meta)
		
		seu_v5@reductions <- seu@reductions
		seu_v5@graphs <- seu@graphs
		seu_v5@neighbors <- seu@neighbors
		seu_v5@misc <- seu@misc
		Idents(seu_v5) <- Idents(seu)
		
	} else if(is(seu$gene, "Assay5")){
		seu_v5 <- seu
	}
	
	# seu_v5 <- seu_v5[,str_detect(seu_v5@meta.data[["phase_level"]], phase_level)]
	
	seu_v5 <- seu_v5[,seu_v5@meta.data[["phase_level"]] == phase_level] 
	
	return(seu_v5)
	
}

merge_seu_paths <- function(seu_paths, phase_level = "g1") {
	seu_paths <- unique(unlist(seu_paths))
	
	seu_paths <- 
		seu_paths |> 
  	set_names(str_extract(seu_paths, "SRR.*(?=_filtered_seu.rds)"))
	
  seu_list <- 
  	seu_paths |> 
  	map(read_seu_by_phase, phase_level = phase_level) |> 
  	identity()
  
  # browser()
  
  marker_genes <- seu_list |>
  	map(~(.x@misc$markers$clusters$presto)) |>
  	map(dplyr::filter, str_detect(Cluster, glue("{phase_level}_[0-9]$"))) |>
  	map(enframe_markers) |>
  	map(slice_head, n = 10) |>
  	map(as_list) |>
  	purrr::list_flatten() |>
  	identity()
  
  # browser()
  
  # features <- SelectIntegrationFeatures(test0)
  # test0 <- PrepSCTIntegration(
  # 	test0,
  # 	anchor.features = features
  # )
  # 
  # # downstream integration steps
  # anchors <- FindIntegrationAnchors(
  # 	test0,
  # 	normalization.method = "SCT",
  # 	anchor.features = features
  # )
  # integrated_object <- IntegrateData(test0, normalization.method = "SCT")
  
  merged_obj<- merge(x=seu_list[[1]], y= unlist(list(seu_list[-1])), merge.dr = FALSE) |>
  	identity()
  
  # browser()
  
  merged_obj <- NormalizeData(merged_obj)
  merged_obj <- FindVariableFeatures(merged_obj)
  merged_obj <- ScaleData(merged_obj)
  merged_obj <- RunPCA(merged_obj)
  
  
  # merged_obj<- merge(x=test0[[1]], y= unlist(list(test0[-1])), merge.dr = TRUE) |>
  # 	identity()
  
  integrated_seu <- IntegrateLayers(
  	object = merged_obj,
  	method = HarmonyIntegration,
  	orig.reduction = "pca",
  	new.reduction = "harmony",
  	verbose = FALSE
  ) |> 
  	JoinLayers() |> 
  	identity()
  # 
  # obj <- IntegrateLayers(
  # 	object = merged_obj,
  # 	method = CCAIntegration,
  # 	orig.reduction = "pca",
  # 	new.reduction = "integrated.cca",
  # 	verbose = FALSE
  # )
  # 
  # obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
  # obj <- FindClusters(obj, resolution = 0.2, cluster.name = "cca_clusters")
  # 
  # obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
  # p1 <- DimPlot(
  # 	obj,
  # 	reduction = "umap.cca",
  # 	group.by = c("sample_id"),
  # 	combine = FALSE, label.size = 2
  # )
  
  return(list("seu" = integrated_seu, "marker_genes" = marker_genes))
  
}

# g1 ------------------------------
# g1_merged_res <- 
	c(
		debranched_seus_1q,
		debranched_seus_2p,
		debranched_seus_6p
		# debranched_seus_16q
	) |> 
	path_file() |> 
	# merge_seu_paths(phase_level = "g1") |> 
	identity()

g1_merged_res$seu$cluster_by_sample <- glue("{g1_merged_res$seu$clusters}_{g1_merged_res$seu$sample_id}")

DimPlot(g1_merged_res$seu, group.by = "cluster_by_sample", reduction = "harmony")

embeddings <-
	g1_merged_res$seu@reductions$harmony@cell.embeddings |>
	data.frame() |>
	tibble::rownames_to_column("cell") |>
	identity()

harmony_centroids <-
	g1_merged_res$seu@meta.data |>
	tibble::rownames_to_column("cell") |>
	dplyr::select(cell, cluster_by_sample) |>
	dplyr::right_join(embeddings) |>
	dplyr::group_by(cluster_by_sample) |>
	summarise(across(starts_with("harmony_"), mean, .names = "mean_{.col}")) |>
	tibble::column_to_rownames("cluster_by_sample") |>
	identity()

cent.tree <- hclust(dist(harmony_centroids), "ward.D2")

plot(cent.tree)
# 
Idents(g1_merged_res$seu) <- g1_merged_res$seu$cluster_by_sample

g1_merged_res$seu$cluster_by_sample <- factor(g1_merged_res$seu$cluster_by_sample, levels = labels(cent.tree))

saveRDS(g1_merged_res, "output/seurat/g1_merged_res.rds")

g1_merged_res <- readRDS("output/seurat/g1_merged_res.rds")

# # s ------------------------------
# s_merged_obj <- merge_seu_paths(
# 	list(
# 		debranched_seus_1q,
# 		debranched_seus_2p,
# 		debranched_seus_6p
# 		# debranched_seus_16q
# 	),
# 	phase_level = "s"
# )
# 
# saveRDS(s_merged_obj, "output/seurat/s_merged_seu.rds")
# s_merged_obj <- readRDS("output/seurat/s_merged_seu.rds")
# 
# # g2 ------------------------------
# g2_merged_obj <- merge_seu_paths(
# 	list(
# 		debranched_seus_1q,
# 		debranched_seus_2p,
# 		debranched_seus_6p
# 		# debranched_seus_16q
# 	),
# 	phase_level = "g2"
# )
# 
# saveRDS(g2_merged_obj, "output/seurat/s_merged_seu.rds")
# g2_merged_obj <- readRDS("output/seurat/g2_merged_seu.rds")


# heatmap_input <- VariableFeatures(g1_merged_res$seu)[1:50]
# heatmap_input <- unlist(g1_merged_res$marker_genes)

# debug(find_all_markers)
seu_scna_col <- g1_merged_res$seu$scna

seu_scna_col[seu_scna_col == ""] <- "diploid"

g1_merged_res$seu$scna <- seu_scna_col

# not aggregated ------------------------------
g1_cluster_plots <- pdf("results/g1_cluster_plots.pdf", height = 12, width = 10)

DimPlot(
	RunUMAP(g1_merged_res$seu, dims = 1:30), group.by = "sample_id", reduction = "umap") +
	labs(title = "uncorrected umap")

heatmap_input <- 
	g1_merged_res$marker_genes |> 
	map(~(.x[1:5])) |> 
	unlist() |>
	identity()

# harmony ------------------------------

DimPlot(g1_merged_res$seu, group.by = "cluster_by_sample", reduction = "harmony")

g1_merged_res$seu <- RunUMAP(g1_merged_res$seu, dims = 1:30)

DimPlot(g1_merged_res$seu, group.by = "cluster_by_sample", reduction = "umap")

# ward_arranged_plot <- 
ggplotify::as.ggplot(seu_complex_heatmap(g1_merged_res$seu,
																				 features = heatmap_input,
																				 group.by = c("phase_level", "cluster_by_sample", "scna"),
																				 # col_arrangement = c("phase_level", "cluster_by_sample", "scna"),
																				 col_arrangement = "ward.D2",
																				 cluster_rows = TRUE,
																				 embedding = "harmony")) +
	labs(title = "not aggregated; ward; harmony")

# pca ------------------------------

DimPlot(g1_merged_res$seu, group.by = "cluster_by_sample", reduction = "pca")

ggplotify::as.ggplot(seu_complex_heatmap(g1_merged_res$seu,
																				 features = heatmap_input,
																				 group.by = c("phase_level", "cluster_by_sample", "scna"),
																				 # col_arrangement = c("phase_level", "cluster_by_sample", "scna"),
																				 col_arrangement = "ward.D2",
																				 cluster_rows = TRUE,
																				 embedding = "pca")) +
	labs(title = "not aggregated; ward; pca")

ggplotify::as.ggplot(seu_complex_heatmap(g1_merged_res$seu,
																				 features = heatmap_input,
																				 group.by = c("phase_level", "cluster_by_sample", "scna"),
																				 col_arrangement = c("phase_level", "cluster_by_sample", "scna"),
																				 # col_arrangement = "ward.D2",
																				 cluster_rows = TRUE,
																				 embedding = "harmony")) +
	labs(title = "not aggregated; cluster_by_sample")

# column_dend(ward_arranged_plot) |> 
# 	labels()

test2 <- find_all_markers(g1_merged_res$seu, metavar = "cluster_by_sample", seurat_assay = "SCT")



# aggregated -----------------------------

test1 <- AggregateExpression(g1_merged_res$seu, group.by = c("ident"), return.seurat = TRUE) |> 
	FindVariableFeatures()

test1$tumor_id <-
	str_extract(test1$orig.ident, "SRR[0-9]*")



# heatmap_input <- VariableFeatures(test1)[1:75]
heatmap_input <- 
	g1_merged_res$marker_genes |> 
	map(~(.x[1:3])) |> 
	unlist() |>
	identity()

ward_arranged_plot <-
	seu_complex_heatmap(test1,
																					 features = heatmap_input,
																					 group.by = c("orig.ident"),
																					 # group.by = c("cluster_by_sample", "phase_level", "scna"),
																					 # col_arrangement = c("cluster_by_sample", "phase_level", "scna"),
																					 col_arrangement = "ward.D2",
																					 cluster_rows = TRUE,
																					 embedding = "harmony")

text_label <- 
column_dend(ward_arranged_plot) |> 
	labels()

test1$orig.ident <- factor(test1$orig.ident, levels = text_label)

ggplotify::as.ggplot(seu_complex_heatmap(test1,
																				 features = heatmap_input,
																				 group.by = c("orig.ident", "tumor_id"),
																				 # group.by = c("cluster_by_sample", "phase_level", "scna"),
																				 # col_arrangement = c("cluster_by_sample", "phase_level", "scna"),
																				 # col_arrangement = c("orig.ident", "tumor_id"),
																				 col_arrangement = c("ward.D2"),
																				 cluster_rows = TRUE,
																				 embedding = "harmony")) + 
	labs(title = "aggregated") + 
	scale_x_discrete(text_label)

test1 <- seurat_reduce_dimensions(test1)

DimPlot(test1, group.by = "orig.ident", label = TRUE)


dev.off()

browseURL("results/g1_cluster_plots.pdf")

g1_merged_res$seu <- find_all_markers(g1_merged_res$seu, metavar = "cluster_by_sample", seurat_assay = "SCT")

heatmap_input <- g1_merged_res$seu@misc$markers$cluster_by_sample$presto |> 
	dplyr::filter(!stringr::str_detect(Gene.Name, "^MT-")) |> 
	dplyr::filter(!stringr::str_detect(Gene.Name, "^RPS")) |> 
	dplyr::filter(!stringr::str_detect(Gene.Name, "^RPL")) |> 
	enframe_markers() |>
	slice_head(n = 5) |>
	unlist() |>
	identity()

Idents(g1_merged_res$seu) <- g1_merged_res$seu$cluster_by_sample

# ward_arranged_plot <- 
	ggplotify::as.ggplot(seu_complex_heatmap(g1_merged_res$seu,
																					 features = heatmap_input,
																					 group.by = c("phase_level", "cluster_by_sample", "scna"),
																					 col_arrangement = c("phase_level", "cluster_by_sample", "scna"),
																					 # col_arrangement = "ward.D2",
																					 cluster_rows = TRUE))

column_dend(ward_arranged_plot) |> 
	labels()





meta_arranged_plot <- ggplotify::as.ggplot(
		seu_complex_heatmap(g1_merged_res$seu,
												features = heatmap_input,
												group.by = c("cluster_by_sample", "phase_level", "scna"),
												col_arrangement = c("cluster_by_sample", "phase_level", "scna"),
												# col_arrangement = "ward.D2",
												cluster_rows = TRUE))
