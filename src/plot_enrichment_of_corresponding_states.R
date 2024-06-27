source("packages.R")
source("functions.R")
library(targets)
library(patchwork)

library('tidyverse')
library('fs')
library('readxl')

tar_load(debranched_seus_2p)

tar_load(debranched_seus_6p)

plot_enrichment_of_corresponding_states <- function(seu_path, clone_comparison = c(4,5)){
	browser()
	seu <- readRDS(seu_path)
	
	seu <- seu[,seu$clone_opt %in% clone_comparison]
	
	enrichment_by_phase <- function(phase_level, seu){
		# browser()
		
		idents <- rev(sort(unique(seu$clone_opt)))
		markers <- 
			seu[,seu$phase_level == phase_level] |> 
			FindMarkers(group.by = "clone_opt", ident.1 = idents[[1]], ident.2 = idents[[2]]) |> 
			enrichment_analysis(gene_set = "hallmark") |> 
			identity()
		
		return(markers)
	}
	
	phase_levels <- janitor::tabyl(seu@meta.data, phase_level, clone_opt) |> 
		identity() |> 
		dplyr::rowwise() |> 
		dplyr::filter(if_all(-phase_level, ~ . > 5)) |> 
		dplyr::pull(phase_level) |> 
		set_names()
	
	map(phase_levels, enrichment_by_phase, seu) |> 
		merge_result() |>
		plot_enrichment() +
		labs(title = seu_path)
	
}

pdf("results/2p_corresponding_terms.pdf", height = 12, width = 12)
# 2p+ ------------------------------

message(debranched_seus_2p[[1]])
plot_enrichment_of_corresponding_states(debranched_seus_2p[[1]], clone_comparison = c(4,5))

message(debranched_seus_2p[[2]])
plot_enrichment_of_corresponding_states(debranched_seus_2p[[2]], clone_comparison = c(2,6))

message(debranched_seus_2p[[3]])
plot_enrichment_of_corresponding_states(debranched_seus_2p[[3]], clone_comparison = c(2,3))

message(debranched_seus_2p[[6]])
plot_enrichment_of_corresponding_states(debranched_seus_2p[[6]], clone_comparison = c(3,4))
dev.off()

browseURL("results/2p_corresponding_terms.pdf")

# 6p+ ------------------------------
pdf("results/6p_corresponding_terms.pdf", height = 12, width = 12)
message(debranched_seus_6p[[1]])
plot_enrichment_of_corresponding_states(debranched_seus_6p[[1]], clone_comparison = c(1,2))

message(debranched_seus_6p[[3]])
plot_enrichment_of_corresponding_states(debranched_seus_6p[[3]], clone_comparison = c(2,3))

dev.off()

browseURL("results/6p_corresponding_terms.pdf")

# plots ------------------------------

# g1

g1_2p_plot <- 
	g1_2p_markers |> 
	plot_enrichment() +
	labs(title= "g1")

# s_g2

s_g2_2p_plot <- 
	s_g2_2p_markers |> 
		plot_enrichment() +
	labs(title= "s_g2")
	
# g2_m
		
		g2_m_2p_plot	 <- 
			g2_m_2p_markers |> 
		plot_enrichment()  +
	labs(title= "g2_m")
	
	# combined plot ------------------------------

merge_result(
	list(
		"g1" = g1_2p_markers,
		"s_g2" = s_g2_2p_markers,
		"g2_m" = g2_m_2p_markers)) |> 
			plot_enrichment()


mrk_tbl <- g1_seu@misc$markers$SCT_snn_res.0.6$presto |> 
	dplyr::filter(!str_detect(Gene.Name, "^MT-")) |> 
	dplyr::filter(!str_detect(Gene.Name, "^RP")) |> 
	dplyr::filter(Adjusted.pvalue < 0.05) |> 
	identity()

mrk_desc <- mrk_tbl |> 
	dplyr::left_join(annotables::grch38, by = c("Gene.Name" = "symbol")) |> 
	identity()

frame_markers <- mrk_tbl |> 
	enframe_markers() |> 
	identity()
	
# asdf ------------------------------

# find_all_markers(g1_seu)

plot_markers(seu, "SCT_snn_res.0.2") +
	plot_markers(g1_seu, "SCT_snn_res.0.2")

g1_seu <- find_all_markers(g1_seu, "SCT_snn_res.0.6", seurat_assay = "SCT", 
													 p_val_cutoff = 1)