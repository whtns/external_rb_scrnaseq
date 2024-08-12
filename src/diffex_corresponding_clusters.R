
library(targets)
source("packages.R")
source("functions.R")

diffex_bw_cluster_set <- function(file_name, w_scna, wo_scna){
	
	browser()
	
	seu <- readRDS(fs::path("output/seurat", file_name))
	
	tumor_id <- str_extract(file_name, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(file_name), "_filtered_seu.*")
	
	w_scna = str_split_1(w_scna, "-")
	wo_scna = str_split_1(wo_scna, "-")
	
	test0 <-
		seu@meta.data |> 
		tibble::rownames_to_column("cell") |> 
		dplyr::mutate(scna_status = dplyr::case_when(
			clusters %in% w_scna ~ "w_scna",
			clusters %in% wo_scna ~ "wo_scna"))
	
	seu$scna_status <- test0$scna_status
	
	seu <-
		seu[,!is.na(seu$scna_status)]
	
	print(sample_id)
	diffex <- FindMarkers(seu, group.by = "scna_status", ident.1 = "w_scna", ident.2 = "wo_scna")
	
	myenrich <- enrichment_analysis(diffex)
	
	enrich_plot <- plot_enrichment(myenrich) + 
		labs(title = sample_id, subtitle = glue("{sample_id} {paste(w_scna, collapse = '-')} v. {paste(wo_scna, collapse = '-')}"))
	
	# newsec ------------------------------
	
	diffex <-
		diffex %>%
		dplyr::mutate(chr = case_when(chr == "X" ~ "23",
																	chr == "Y" ~ "24",
																	TRUE ~ as.character(chr))) %>%
		dplyr::mutate(chr = str_pad(chr, side = "left", pad = "0", width = 2)) %>%
		dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", " ")) %>%
		tibble::rownames_to_column("symbol") %>%
		dplyr::mutate(rownames = symbol) %>%
		tibble::column_to_rownames("rownames")
	
	mytitle = sample_id
	mysubtitle = mysubtitle
	
	ref_var <-
		diffex$chr %>%
		set_names(.)
	
	chrs <- str_pad(as.character(1:24), side = "left", pad = "0", width = 2)
	
	mypal <- scales::hue_pal()(length(chrs))
	names(mypal) <- chrs
	
	custom_cols <- mypal[ref_var]
	
	FCcutoff = summary(abs(diffex$avg_log2FC))[[5]]
	
	selected_genes <-
		diffex %>%
		dplyr::filter(abs(avg_log2FC) > 0.05, p_val_adj < 0.1) %>%
		dplyr::pull(symbol)
	
	
	volcano_plot <- EnhancedVolcano(diffex,
														lab = rownames(diffex),
														selectLab = selected_genes,
														labSize = 4,
														x = 'avg_log2FC',
														y = 'p_val_adj',
														FCcutoff = FCcutoff,
														pCutoff = 5e-2,
														colCustom = custom_cols,
														max.overlaps = 25) +
		aes(color = chr) +
		# facet_wrap(~chr) +
		labs(title = mytitle, subtitle = mysubtitle)
	
	# endsec ------------------------------
	
	return(list("enrich" = enrich_plot, "volcano" = volcano_plot))
	
}

corresponding_states_dictionary = 
					 tibble::tribble(
                                          ~file_name, ~w_scna, ~wo_scna,
        "SRR13884246_branch_5_filtered_seu_2p.rds", "g1_1-g1_4-g1_0",      "g1_2",
        "SRR13884246_branch_5_filtered_seu_2p.rds",     "s_g2_7",      "s_6",
           "SRR13884247_branch_6_filtered_seu.rds", "g1_0-g1_2-g1_5",    "g1_3-g1_6",
                 "SRR13884248_filtered_seu_2p.rds",   "g1_0-g1_1",    "g1_6-g1_7",
                 "SRR17960484_filtered_seu_2p.rds",     "g1_1",    "g1_3-g1_4",
                 "SRR17960484_filtered_seu_2p.rds",     "s_g2_5",      "s_g2_0",
                    "SRR13884247_filtered_seu.rds",   "g1_0-g1_1",    "g1_4-g1_5",
                 "SRR17960484_filtered_seu_6p.rds",     "g1_2",      "g1_3"
        )


test0 <- pmap(corresponding_states_dictionary, diffex_bw_cluster_set)


pdf("results/corresponding_clusters.pdf")
print(test0)
dev.off()

browseURL("results/corresponding_clusters.pdf")
