#!/usr/bin/env Rscript

source("packages.R")
source("functions.R")

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)

seus <- c(
	"output/seurat/SRR14800534_filtered_seu.rds",
	"output/seurat/SRR14800535_filtered_seu.rds",
	"output/seurat/SRR14800536_filtered_seu.rds"
) |> 
	set_names()

seus <- set_names(seus, str_extract(names(seus), "SRR[0-9]*"))

seus <- seus |> 
	map(readRDS)


subset_seu_by_clones <- function(seu, clones){
	seu <- seu[,seu$clone_opt %in% clones]
	
	return(seu)
}

seus <- map2(seus, list(
	c(1,2),
	c(1,2),
	c(1,2)
	),
	subset_seu_by_clones)

integrated_seu <- seuratTools::integration_workflow(seus)

saveRDS(integrated_seu, "output/seurat/integrated_seu_16q_trio.rds")

integrated_seu <- readRDS("output/seurat/integrated_seu_16q_trio.rds")

integrated_seu <- ScaleData(integrated_seu)

integrated_seu <- seurat_reduce_dimensions(integrated_seu)

integrated_seu <- seurat_cluster(integrated_seu, resolution = seq(0.2, 2.0, by = 0.2), seurat_assay = "integrated")

saveRDS(integrated_seu, "output/seurat/integrated_seu_16q_complete.rds")


make_cc_plot(integrated_seu, var_y = "integrated_snn_res.0.4")




