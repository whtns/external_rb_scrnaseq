#!/usr/bin/env Rscript

library(singlet)
library(Seurat)
library(dplyr)
library(ggplot2)
library(fs)
library(tidyverse)
library(getopt)
library(sceasy)
library(glue)
library(fs)
library(purrr)
library(tidyverse)
library(seuratTools)
library(fastcluster)

# seu_path = "output/seurat/SRR14800534_filtered_seu.rds"
# pdf_path = "output/mosaicmpi/SRR14800534_heatmaps.pdf"

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}

sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")

tumor_id <- str_extract(seu_path, "SRR[0-9]*")

seu <- readRDS(seu_path)

usages = dir_ls(glue("output/mosaicmpi/{sample_id}/"), glob = "*usage_k*.txt") %>% 
	set_names(str_extract_all(., "(?<=k)[0-9]*")) %>% 
	map(read_tsv) %>% 
	map(tibble::column_to_rownames, "...1") %>% 
	imap(~set_names(.x, paste0(.y, "_", colnames(.x)))) %>% 
	identity()

usages <-
	usages[as.character(sort(as.numeric(names(usages))))]

colnames(usages[["6"]]) <- paste0('factor.', 1:6)

seu0 <- AddMetaData(seu, bind_cols(usages))

seu0$sample <- sample_id

# programs <- read_tsv(glue("output/mosaicmpi/{sample_id}/_programs.csv")) %>% 
# 	dplyr::rename(symbol = `...1`) %>% 
# 	tidyr::pivot_longer(-symbol, names_to = "factor", values_to = "usage") %>% 
# 	dplyr::filter(str_detect(factor, "^6")) %>% 
# 	dplyr::rename(Gene.Name = symbol, Cluster = factor) %>% 
# 	dplyr::group_by(Cluster) %>% 
# 	dplyr::arrange(desc(usage)) %>% 
# 	dplyr::slice_head(n = 6) %>% 
# 	identity()

pdf(pdf_path, height = 10, width = 10)

if(all(c("clusters", "scna") %in% colnames(seu0@meta.data))){
	
	# factor_plot <- seu_gene_heatmap(seu0, marker_col = "clusters", 
	# 																group.by = c("clusters", "scna", "S.Score", "G2M.Score", colnames(usages[["6"]])), 
	# 																col_arrangement = c("clusters", "scna"), 
	# 																column_split = "clusters", 
	# 																hide_legends = colnames(usages[["6"]]),
	# 																heatmap_features = programs) + labs(title = sample_id, subtitle = "6")
	# print(factor_plot)
	
	cluster_plot <- seu_gene_heatmap(seu0, marker_col = "clusters", 
																	 group.by = c("clusters", "scna", "S.Score", "G2M.Score", colnames(usages[["6"]])), 
																	 col_arrangement = c("clusters", "scna"), 
																	 column_split = "clusters", 
																	 hide_legends = colnames(usages[["6"]])) + labs(title = sample_id, subtitle = "6")
	print(cluster_plot)
	
	
	
	# cluster_plots <- map(names(usages), ~seu_factor_heatmap(seu0, group.by = c("clusters", "scna", "sample"), col_arrangement = c("clusters", "scna"), factor_cols = colnames(usages[[.x]]), column_split = "clusters")[["plot"]] + labs(title = sample_id, subtitle = .x))
	# print(cluster_plots)
	
	# arranged_plots <- map(names(usages), ~seu_factor_heatmap(seu0, group.by = c("clusters", "scna"), col_arrangement = "ward.D2", factor_cols = colnames(usages[[.x]]))[["plot"]] + labs(title = sample_id, subtitle = .x))
	# print(arranged_plots)
}  else {
	cluster_plots <- map(names(usages), ~seu_factor_heatmap(seu0, group.by = c("SCT_snn_res.0.6"), col_arrangement = c("SCT_snn_res.0.6"), factor_cols = colnames(usages[[.x]]), column_split = "SCT_snn_res.0.6")[["plot"]] + labs(title = sample_id, subtitle = .x))
	print(cluster_plots)
	# arranged_plots <- map(names(usages), ~seu_factor_heatmap(seu0, group.by = c("SCT_snn_res.0.6", "scna"), col_arrangement = "ward.D2", factor_cols = colnames(usages[[.x]]))[["plot"]] + labs(title = sample_id, subtitle = .x))
	# print(arranged_plots)
}

dev.off()
