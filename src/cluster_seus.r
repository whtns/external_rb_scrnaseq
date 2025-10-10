source("packages.R"); source("functions.R")

seu_paths <- c(
	"output/seurat/SRR27187900_seu.rds",
	"output/seurat/SRR27187901_seu.rds"
)

# undebug(find_all_markers)

run_clustering <- function(seu_path) {
	seu <- readRDS(seu_path)
	DefaultAssay(seu) <- "gene"
	seu <- seuratTools::clustering_workflow(seu, resolution = c(0.2, 0.4))
	DefaultAssay(seu) <- "SCT"
	saveRDS(seu, seu_path)
	
	cluster_nums <- seu@misc$markers$gene_snn_res.0.2$presto |> 
			dplyr::mutate(Cluster = as.numeric(Cluster)) |> 
			dplyr::group_by(Cluster) |> 
			dplyr::slice_head(n = 5) |> 
			dplyr::ungroup()  |> 
			with(table(Cluster))
	
	return(cluster_nums)
}

run_clustering(seu_paths[2])

cluster_nums <- map(seu_paths, run_clustering)

seu@misc$markers$gene_snn_res.0.2$presto |> 
			dplyr::mutate(Cluster = as.numeric(Cluster)) |> 
			dplyr::group_by(Cluster) |> 
			dplyr::slice_head(n = 5) |> 
			dplyr::ungroup()  |> 
			with(table(Cluster))


View(find_all_markers)

