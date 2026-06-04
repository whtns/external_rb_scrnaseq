#!/usr/bin/env Rscript

library(targets)
library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)

source("packages.R")
source("functions.R")

make_clone_plot <- function(seu, var_x = "clusters", var_y = "scna"){
	# browser()
	
	seu_meta <-
		seu@meta.data
	
	summarized_clones <-
		seu_meta %>%
		dplyr::select(.data[[var_x]], .data[[var_y]]) %>%
		dplyr::mutate({{var_x}} := "all")
	
	clone_plot <- 
		seu_meta |> 
		tibble::rownames_to_column("cell") |> 
		ggplot() + 
		geom_bar(data = summarized_clones, position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]])) +
		geom_bar(position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]])) +
		scale_x_discrete(limits = c(rev(levels(seu_meta[[var_x]])), "all")) +
		coord_flip() +
		facet_wrap(~sample_id, scales = "free_y") + 
		NULL
	
	cc_data <- FetchData(seu, c(var_y, var_x, "G2M.Score", "S.Score", "Phase", "scna"))
	
	centroid_data <-
		cc_data %>%
		dplyr::group_by(.data[[var_x]]) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		dplyr::mutate({{var_x}} := factor(.data[[var_x]], levels = levels(cc_data[[var_x]]))) %>%
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	centroid_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[var_x]], color = .data[[var_x]])) +
		geom_point(size = 0.1) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[var_x]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
		guides(fill = "none", color = "none") +
		NULL
	
	facet_cell_cycle_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[var_x]], color = .data[[var_y]])) +
		geom_point(size = 0.1) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[var_x]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
		facet_wrap(~ .data[[var_x]], ncol = 2) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		NULL
	
	# mypatch <- wrap_plots(list(clone_plot, facet_cell_cycle_plot)) + 
	# 	plot_layout(ncol = 1)
	
	return(list(clone_plot, facet_cell_cycle_plot))
}

clone_dist_before_and_after <- function(seu, resolution = 0.2){
	
	layout =
	"
	AACC
	AACC
	BBDD
	BBDD
	"
	
	test0 <-	list(
		make_clone_plot(seus[[1]], var_x = glue("old_SCT_snn_res.{resolution}")),
		make_clone_plot(seus[[1]], var_x = glue("new_SCT_snn_res.{resolution}"))) |> 
		purrr::flatten() |> 
		wrap_plots() +
		plot_layout(design = layout)
	
}

tar_load(c("debranched_seus_2p", "debranched_seus_6p", "large_clone_comparisons"))

debug(seurat_integration_pipeline)

retrieve_clones_of_interest <- function(seu, sample_id, large_clone_comparisons, scna_of_interest){
	clone_comparisons <- names(large_clone_comparisons[[sample_id]])
	clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
	retained_clones <- clone_comparison %>%
		str_extract("[0-9]_v_[0-9]") %>%
		str_split("_v_", simplify = TRUE)
	
	seu <- seu[,seu$clone_opt %in% retained_clones]
	
	return(seu)
	
}

stash_colnames <- function(seu, scna_of_interest = "6p\\+"){
	sct_colnames <- str_subset(colnames(seu@meta.data), "^SCT")
	
	old_colnames <- str_replace(sct_colnames, "SCT", "old_SCT")
	new_colnames <- str_replace(sct_colnames, "SCT", "new_SCT")
	
	seu@meta.data[old_colnames] <- seu@meta.data[sct_colnames]
	
	seu <- seurat_cluster(seu, resolution = seq(0.2, 2.0, by = 0.2))
	
	seu@meta.data[new_colnames] <- seu@meta.data[sct_colnames]
	
	seu@meta.data[sct_colnames] <- seu@meta.data[old_colnames]
	
	seu$scna <- factor(ifelse(str_detect(seu$scna, scna_of_interest), "w_scna", "wo_scna"), levels = c("wo_scna" , "w_scna"))
	
	for(colname in new_colnames){
		seu <- find_all_markers(seu, metavar = colname)
	}
	
	for(colname in old_colnames){
		seu <- find_all_markers(seu, metavar = colname)
	}
	
	return(seu)
}

seus <- map(debranched_seus_6p, readRDS) |> 
	imap(retrieve_clones_of_interest, large_clone_comparisons, scna_of_interest = "6p\\+") |> 
	map(stash_colnames) |> 
	# integration_workflow(find_markers = FALSE) |> 
	# saveRDS("output/seurat/integrated_6p/seurat_6p_integrated.rds") |> 
	identity()

map2(seus, debranched_seus_6p, saveRDS)


test1 <- 
	seq(0.2, 1.0, by = 0.2) |> 
	set_names() |> 
	map(~clone_dist_before_and_after(seu, .x))

pdf("results/")



# 2p ------------------------------

seus <- map(debranched_seus_2p, readRDS) |> 
	imap(retrieve_clones_of_interest, large_clone_comparisons, scna_of_interest = "2p") |> 
	integration_workflow(find_markers = FALSE) |> 
	saveRDS("output/seurat/integrated_2p/seurat_2p_integrated.rds")


# 6p ------------------------------
seus <- map(debranched_seus_6p, readRDS) |> 
	imap(retrieve_clones_of_interest, large_clone_comparisons, scna_of_interest = "2p") |> 
	integration_workflow() |> 
	saveRDS("output/seurat/integrated_6p/seurat_6p_integrated.rds")
