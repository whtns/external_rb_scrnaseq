#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
source("packages.R")
source("functions.R")


pull_scna_segments <- function(nb_path = "output/numbat_sridhar/SRR14800534_numbat.rds", chrom = "1"){
	mynb <- readRDS(nb_path)
	
	segments <- mynb$clone_post %>%
		dplyr::left_join(mynb$joint_post, by = "cell") %>%
		# dplyr::filter(clone_opt %in% idents) %>%
		# dplyr::filter(seg %in% mysegs) %>%
		dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
		dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
		dplyr::filter(!cnv_state_map == "neu") %>%
		dplyr::filter(seqnames == chrom) |> 
		plyranges::as_granges() %>%
		identity()
}

diffex_per_cluster <- function(seu, ident_1 = "w_scna", ident_2 = "wo_scna", mycluster = 'g1_1', segments = NULL){
	seu <- seu[,seu$clusters %in% mycluster]
	
	diffex <- FindMarkers(seu, ident.1 = ident_1, ident.2 = ident_2, group = "scna", test.use = "wilcox", assay = "gene", logfc.threshold = 0.1) |> 
		tibble::rownames_to_column("symbol") %>%
		dplyr::left_join(annotables::grch38, by = "symbol") %>%
		dplyr::distinct(ensgene, .keep_all = TRUE) %>%
		dplyr::mutate(seqnames = chr) %>%
		dplyr::filter(!is.na(start), !is.na(end)) %>%
		plyranges::as_granges() |>
		identity()
	
	# cis 
	
	cis_diffex <- 
		diffex |> 
		plyranges::join_overlap_intersect(segments) %>%
		as_tibble() %>%
		# dplyr::mutate(log2_sign = dplyr::case_when(
		# 	cnv_state_map == "amp" ~ -1,
		# 	cnv_state_map == "del" ~ 1
		# )) %>%
		# dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
		dplyr::select(-any_of(c(
			"CHROM", "seg_start",
			"seg_end", "cnv_state_map", "log2_sign"
		))) %>%
		dplyr::filter(!str_detect(chr, "CHR_")) %>%
		dplyr::distinct(symbol, .keep_all = TRUE)
	
	return(cis_diffex)
	
}



# 16q  ------------------------------

segs_16q <- sapply(c("output/numbat_sridhar/SRR14800534_numbat.rds",
										 "output/numbat_sridhar/SRR14800535_numbat.rds",
										 "output/numbat_sridhar/SRR14800536_numbat.rds"), pull_scna_segments, chrom = "16") |> 
	as("GRangesList") |> 
	unlist() |> 
	reduce_ranges()

seu <- readRDS("output/seurat/integrated_16q/integrated_seu_16q_complete.rds")

myclusters <- seu$clusters |> 
	levels() |> 
	set_names(identity)

diffex_16q <- map(myclusters, ~diffex_per_cluster(seu, ident_1 = "16q-", ident_2 = "", mycluster = .x, segments = segs_16q)) |> 
	dplyr::bind_rows(.id = "cluster")

filtered_diffex_16q <- 
	diffex_16q |> 
	dplyr::filter(p_val_adj <= 0.1, seqnames == 16) |> 
	dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
	# dplyr::filter(avg_log2FC > 0.1) |> 
	dplyr::filter(!cluster %in% c("hsp_8", "hypoxia_2")) |>
	dplyr::arrange(desc(cluster), desc(avg_log2FC)) |> 
	dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |> 
	identity()

diffex_16q |> 
	dplyr::filter(p_val_adj <= 0.1) |> 
	janitor::tabyl(cluster)

ggplot(filtered_diffex_16q, aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) + 
	geom_point() + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
	scale_color_gradient2() +
	labs(title = "integrated 16q- diffex per cluster")

ggsave("results/integrated_16q_diffex_by_cluster.pdf", w = 4, h = 4) |> 
	browseURL()
