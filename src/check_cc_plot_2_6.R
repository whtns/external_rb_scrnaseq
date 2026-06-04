#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(Seurat)
source("packages.R")
source("functions.R")

make_clone_plot <- function(seu, var_y = "clusters", var_x = "scna"){
	# browser()
	clone_plot <- 
		seu@meta.data |> 
		tibble::rownames_to_column("cell") |> 
		ggplot() + 
		geom_bar(position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
		# geom_bar(data = summarized_clusters, position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
		scale_x_discrete(limits = rev) +
		coord_flip() +
		facet_wrap(~sample_id, scales = "free_y") + 
		NULL
	
	return(clone_plot)
}

scna_seu <- readRDS("output/seurat/integrated_2p/seurat_2p_integrated.rds")
scna_of_interest = "2p\\+"

scna_seu@meta.data[,str_detect(colnames(scna_seu@meta.data), "sample_id")] <- NULL

scna_seu$sample_id <- scna_seu$batch

scna_seu$old_scna <- scna_seu$scna 

scna_seu$scna <- ifelse(str_detect(scna_seu$old_scna, scna_of_interest), "w_scna", "wo_scna")


test0 <- map(glue("integrated_snn_res.{seq(0.2, 2.0, by = 0.2)}"), ~make_cc_plot(scna_seu, var_y = .x))


test2 <- map(glue("integrated_snn_res.{seq(0.2, 2.0, by = 0.2)}"), ~make_clone_plot(scna_seu, var_x = "scna", var_y = .x) + 
						 	facet_wrap(~sample_id, nrow = 1) + 
						 	labs(title = .x))

plot_path <- glue("results/{scna_of_interest}_integrated_ccplots.pdf")

pdf(plot_path, w = 12, h = 8)
		map(test0, 2)
		print(test2)
dev.off()

browseURL(plot_path)



# diffex <- FindMarkers(scna_seu, group.by = "scna", ident.1 = "w_scna", ident.2 = "wo_scna") |> 
# 	tibble::rownames_to_column("symbol") |> 
# 	dplyr::filter(p_val_adj < 0.05) |> 
# 	dplyr::arrange(desc(abs(avg_log2FC))) |> 
# 	identity()
# 
# VlnPlot(scna_seu, features = diffex$symbol[1:6], group.by = "scna")


# DimPlot(scna_seu, group.by = "scna")

# saveRDS(scna_seu, "output/seurat/integrated_2p/seurat_2p_integrated.rds")