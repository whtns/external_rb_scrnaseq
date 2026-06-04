library(tidyverse)
library(fs)
library(seuratTools)
library(targets)
library(patchwork)

tar_load(c("debranched_seus_16q"))

myseus <-
	debranched_seus_16q %>% 
	set_names(path_file(str_remove(., "_filtered_seu.rds"))) %>%
	identity()


myviolins <- map2(myseus, large_clone_comparisons[names(myseus)], ~{
	# browser()
	
	sample_id <-path_file(str_remove(.x, "_filtered_seu.rds"))
	
	# print(sample_id)
	seu <- readRDS(.x)
	
	# clone_comparison <- names(.y) %>% 
	# 	tibble::enframe("row_num", "clones") %>%
	# 	dplyr::filter(str_detect(clones, ".*16q.*")) %>% 
	# 	dplyr::mutate(clones = str_extract(clones, "[0-9]_v_[0-9]")) %>%
	# 	dplyr::mutate(clones = str_split(clones, "_v_")) %>%
	# 	tidyr::unnest(clones) %>% 
	# 	dplyr::pull(clones) %>% 
	# 	identity()
	# 
	# # print(clone_comparison)
	# 
	# seu <- seu[,seu$clone_opt %in% clone_comparison]
	
	# RBL2_seu <- seu[,FetchData(seu, "RBL2") > 0]
	# SETD6_seu <- seu[,FetchData(seu, "SETD6") > 0]
	
	mean_col_w_points <- function(seu, plot_var = "scna", features= "RBL2", split.by = NULL){
		# browser()
		DefaultAssay(seu) <- "gene"
		
		test0 <- FetchData(seu, c("scna", "clusters", features), layer = "data")
		
		fill_var = plot_var
		
		if(!is.null(split.by)){
			fill_var = split.by
		}
		
		ggplot(test0, aes(x = .data[[plot_var]], y = .data[[features]], fill = .data[[fill_var]], color = .data[[fill_var]], group = .data[[fill_var]])) +
			stat_summary(aes(y = .data[[features]]), fun = "mean", geom = "bar", position = "dodge") + 
			# geom_col(position = "dodge", stat = "summary", fun.y = "mean") +
			geom_point(position = position_dodge(width = .9)) +
			# scale_y_log10() +
			# geom_violin() + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
			NULL
	}
	
	mean_col_w_points(seu, plot_var = "clusters", features = c("RBL2"), split.by = "scna")
	
	# getplot_violin(seu, plot_var = "scna", features = c("RBL2"), log = TRUE) + 
	# 	plot_violin(seu, plot_var = "clusters", features = c("RBL2"), split.by = "scna", log = TRUE) + 
	# 	plot_violin(seu, plot_var = "scna", features = c("SETD6"), log = TRUE) + 
	# 	plot_violin(seu, plot_var = "clusters", features = c("SETD6"), split.by = "scna", log = TRUE) + 
	# 	plot_layout(ncol = 2, nrow = 2) + 
	# 	plot_annotation(title = sample_id)
	
	mean_col_w_points(seu, plot_var = "scna", features = c("RBL2")) + 
		mean_col_w_points(seu, plot_var = "clusters", features = c("RBL2"), split.by = "scna") + 
		mean_col_w_points(seu, plot_var = "scna", features = c("SETD6")) + 
		mean_col_w_points(seu, plot_var = "clusters", features = c("SETD6"), split.by = "scna") + 
		plot_layout(ncol = 2, nrow = 2) + 
		plot_annotation(title = sample_id)
	
	
})

View(Seurat.FindMarkers)



pdf("results/RBL2_SETD6_variation_across_clones_clusters_in_16qloss_samples.pdf", width = 10, height = 10)
print(myviolins)
dev.off()

browseURL("results/RBL2_SETD6_variation_across_clones_clusters_in_16qloss_samples.pdf")


test0 <- FindMarkers(seu, group.by= "scna", ident.1 = "", ident.2 = "16q-", assay = "gene", min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)
