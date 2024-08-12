#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

seu <- readRDS("output/seurat/SRR17960484_filtered_seu_2p.rds")

source("packages.R")
source("functions.R")

debug(find_cluster_pairwise_distance)
# debug(calculate_avg_cluster_distance)

test0 <- find_cluster_pairwise_distance(seu, "new_SCT_snn_res.0.6") |> 
	dplyr::mutate(similarity = 1-avg_dist) |> 
	dplyr::arrange(x) |> 
	dplyr::mutate(x = fct_inseq(x), y = fct_inseq(y)) |> 
	identity()



ggplot(test0, aes(x = x, y = y, fill = similarity)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = fivenum(test0$similarity)[3], limits = c(0.2,0.8), oob=scales::squish) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(test0) +
	geom_col(aes(x = x, y = similarity, fill = y), position = position_dodge()) + 
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


test0 |>
	arrange(x, desc(similarity)) |>
	mutate(group_order = forcats::fct_inorder(interaction(x, y))) |> 
	ggplot(mapping = aes(x = x, y = similarity, fill = y, group = group_order)) +
	geom_col(position = position_dodge()) + 
	geom_text(aes(label = y),
						position = position_dodge(width = 0.9), angle = 90, hjust = 1
	)
