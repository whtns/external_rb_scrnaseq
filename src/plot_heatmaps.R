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

seu_path = "output/seurat/SRR27187899_filtered_seu.rds"

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

seu0 <- AddMetaData(seu, bind_cols(usages))

source("functions.R")

# debug(seu_factor_heatmap)


# heatmaps <- plot_seu_marker_heatmap(seu0, sample_id = sample_id, tumor_id = tumor_id, 
# 																		othercols = colnames(usages$k90),
# 																		use_raster = TRUE)


pdf(pdf_path, height = 10, width = 10)
map(names(usages), ~ggplotify::as.ggplot(seu_factor_heatmap(seu0, group.by = c("clusters", "scna"), col_arrangement = c("clusters", "scna"), factor_cols = colnames(usages[[.x]]), column_split = TRUE)) + labs(title = sample_id, subtitle = .x))

dev.off()
