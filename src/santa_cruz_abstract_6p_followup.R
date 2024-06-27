source("packages.R")
source("functions.R")
library(targets)

seu <- readRDS("output/seurat/SRR17960484_filtered_seu.rds")

seu$clusters


mydiffex <- Seurat::FindMarkers(seu, group.by = "clusters", ident.1 = "g1_2", ident.2 = "g1_3") %>% 
	tibble::rownames_to_column("symbol") %>% 
	dplyr::left_join(annotables::grch38, by = "symbol") %>% 
	dplyr::select(description, everything()) %>% 
	dplyr::distinct(symbol, .keep_all = TRUE) %>% 
	dplyr::filter(abs(avg_log2FC) > 1) %>% 
	dplyr::filter(p_val_adj < 0.05) %>% 
	dplyr::arrange(avg_log2FC, p_val_adj) %>% 
	identity()

down_csv(mydiffex)
