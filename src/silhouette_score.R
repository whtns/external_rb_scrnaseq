# silhouette metric

source("packages.R")
source("functions.R")
library(targets)

seu <- readRDS("output/seurat/SRR14800534_filtered_seu.rds")

library(cluster, quietly = TRUE)

calc_silhouette <- function(seu_path, assay = "SCT") {
	
	seu <- readRDS(seu_path)
	
	resolutions <- glue("{assay}_snn_res.{seq(0.2, 2.0, by = 0.2)}") %>%
		set_names(.)
	
	reduction = "pca"
	dims = 1:30
	
	for(resolution in resolutions){
		dist.matrix <- dist(x = Embeddings(object = seu[[reduction]])[, dims])
		clusters <- seu@meta.data[[resolution]]
		sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
		
		sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")
		
		seu[[sil_res]] <- sil[, 3]	
	}
	
	for(resolution in resolutions){
		sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")
		
		seu_meta <- 
			seu@meta.data |> 
			tibble::rownames_to_column("cell") |> 
			dplyr::arrange(.data[[resolution]], desc(.data[[sil_res]])) |> 
			dplyr::mutate(cell = factor(cell, levels = unique(cell)))
		
		mean_sil <- mean(seu_meta[[sil_res]])
		
		plot_list[[sil_res]] = ggplot(data = seu_meta, aes(x = .data[["cell"]], y = .data[[sil_res]], fill = .data[[resolution]])) + 
			geom_col(outlier.size = 0.1)  +
			geom_hline(aes(yintercept = mean_sil)) + 
			xlab("Method") + ylab("Silhoutte Metric") +
			labs(title = resolution)
	}
	
	plot_path <- 
  
  return(seu)
}

seu0 <- calc_silhouette(seu, "SCT")



plot_silhouette <- function(seu, assay = "SCT") {
	
	resolutions <- glue("{assay}_snn_res.{seq(0.2, 2.0, by = 0.2)}") %>%
		set_names(.)
	
	plot_list <- list()
	for(resolution in resolutions){
		sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")
	
		seu_meta <- 
			seu@meta.data |> 
			tibble::rownames_to_column("cell") |> 
			dplyr::arrange(.data[[resolution]], desc(.data[[sil_res]])) |> 
			dplyr::mutate(cell = factor(cell, levels = unique(cell)))
		
		mean_sil <- mean(seu_meta[[sil_res]])
		
	  plot_list[[sil_res]] = ggplot(data = seu_meta, aes(x = .data[["cell"]], y = .data[[sil_res]], fill = .data[[resolution]])) + 
	  	geom_col(outlier.size = 0.1)  +
	  	geom_hline(aes(yintercept = mean_sil)) + 
	  	xlab("Method") + ylab("Silhoutte Metric") +
	  	labs(title = resolution)
	}
	
	return(plot_list)
}

plot_silhouette(seu0, "SCT")

test0 <- map(resolutions, ~plot_silhouette(seu0, .x))

pdf("~/tmp/534_silhouettes.pdf")
test0
dev.off()

browseURL("~/tmp/534_silhouettes.pdf")

plot_silhouette(seu, "SCT_snn_res.0.6") + 
	plot_silhouette(seu, "SCT_snn_res.0.4") + 
	plot_silhouette(seu, "SCT_snn_res.0.2") + 
	plot_layout(ncol =1)

# mixing metric
max.k <- 300
mm <- max.k - MixingMetric(object = seu, grouping.var = "replicate", reduction = reduction, dims = dims, max.k = max.k)

DefaultAssay(object = dataset) <- "RNA"
# Local structure preservation
ls <- LocalStruct(object = dataset, grouping.var = "replicate", reduction = reduction, reduced.dims = dims, orig.dims = 1:30)
ls <- unname(obj = unlist(x = ls))

all.metrics <- list(
	silhouette = dataset$sil, 
	mixing.metric = mm,
	local.struct = ls
)

out <- capture.output(FindClusters(seu), type =  "output")
modularity <- as.numeric(unlist(strsplit(x = out[7], split = ":"))[2])
