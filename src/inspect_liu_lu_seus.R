#!/usr/bin/env Rscript

source("packages.R")
source("functions.R")

library('tidyverse')
library('fs')
library('readxl')
library(patchwork)




liu_supp_tables <- read_liu_lu_supp_tables()

liu_meta <- 
	liu_supp_tables$`Fig 5c` %>% 
	dplyr::mutate(sample = str_remove(cell, "[ACGT]*-1_")) %>% 
	tibble::column_to_rownames("cell") %>% 
	identity()


mymat <- Seurat::Read10X("data/GSE249995/SRR27187899/")

liu_lu_seu_paths <- 
	c(
		"output/seurat/SRR27187899_seu.rds",
		"output/seurat/SRR27187900_seu.rds",
		# "output/seurat/SRR27187901_seu.rds",
		"output/seurat/SRR27187902_seu.rds"
	)

for(seu_path in liu_lu_seu_paths){
	print(seu_path)
	seu <- readRDS(seu_path)
	seu$gene_snn_res.0.2 <- seu$SCT_snn_res.0.2
	saveRDS(seu, seu_path)
}

liu_lu_seus <- 
	liu_lu_seu_paths %>% 
	set_names(str_extract(., "SRR[0-9]*")) %>% 
	map(readRDS) %>% 
	map(~{.x[["gene_snn_res.0.2"]]})

test0 <- 
liu_lu_seus %>% 
	dplyr::bind_rows(.id = "sample_id") %>% 
	dplyr::select(sample_id, gene_snn_res.0.2) %>% 
	dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>% 
	dplyr::distinct() %>% 
	dplyr::mutate(abbreviation = as.character(gene_snn_res.0.2), remove = NA) %>% 
	dplyr::arrange(sample_id, abbreviation)



patch_plots <- function(seu, sample_id = "sample_id", metavar = "SCT_snn_res.0.2"){
	plot_cc(seu, group.by = metavar, color.by = metavar) + 
	plot_markers(seu, metavar = metavar, seurat_assay = "SCT") + 
		patchwork::plot_annotation(title = sample_id)
}

seu_plots <- imap(liu_lu_seus, patch_plots)

pdf("results/liu_lu_seu_annotation.pdf", width = 10)
seu_plots
dev.off()

browseURL("results/liu_lu_seu_annotation.pdf")

# test0 ------------------------------

plot_cc(liu_lu_seus[[1]], group.by = "SCT_snn_res.1.6", assay = "SCT")
