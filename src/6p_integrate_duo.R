
library(Seurat)
library(seuratTools)
library(patchwork)
library(glue)

seu_899 <- readRDS("output/seurat/SRR27187899_filtered_seu.rds")
seu_899$batch <- c("SRR27187899")
seu_4248 <- readRDS("output/seurat/integrated_6p/SRR13884248_integrated_6p_filtered_seu.rds")
seu_4248$batch <- c("SRR13884248")

# seu <- readRDS("output/seurat/integrated_6p/integrated_seu_6p_complete.rds")

seu1 <- 
	seurat_integrate(
		list(
			"SRR27187899" = seu_899, 
			"SRR13884248" = seu_4248)
	)

DefaultAssay(seu1) <- "gene"

seu1 <- JoinLayers(seu1)

DefaultAssay(seu1) <- "integrated"

seu2 <- seuratTools::seurat_cluster(seu1, resolution = seq(0.2, 2.0, by = 0.2)) |> 
	find_all_markers(seurat_assay = "integrated", metavar = glue("integrated_snn_res.0.{seq(2,6,2)}"))

plot_markers(seu2, metavar = "integrated_snn_res.0.6")

seu2$scna <- factor(seu2$scna, levels = c("wo_scna", "w_scna"))
	
saveRDS(seu2, "output/seurat/integrated_6p/integrated_seu_6p_duo2.rds")
	
test0 <- make_clone_distribution_figure(
	seu_path = "output/seurat/integrated_6p/integrated_seu_6p_duo2.rds", 
	group.bys = glue("integrated_snn_res.0.{seq(2,6,2)}"))

# set new clusters ------------------------------

tar_load(c("cluster_orders"))

debug(assign_phase_clusters)

seu0 <- 
	assign_phase_clusters(seu = seu1, 
												file_id = "integrated_seu_6p_duo.rds", 
												cluster_orders = cluster_orders)

saveRDS(seu0, "output/seurat/integrated_6p/integrated_seu_6p_duo.rds")



