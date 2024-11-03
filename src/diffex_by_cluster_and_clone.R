library(targets)
source("packages.R")
source("functions.R")

select_recurrent_genes <- function(gene_list, x_var = "sample_id", recurrence_threshold = 2){
	# browser()
	gene_list <- 
		gene_list |> 
		dplyr::filter(recurrence >= recurrence_threshold) |> 
		dplyr::arrange(desc(abs_mean_FC)) |> 
		identity()
	
	gene_list$symbol <- factor(gene_list$symbol, levels = unique(gene_list$symbol))
	
	top_genes <- unique(as.character(gene_list$symbol))[1:50]
	
	gene_list |> 
		dplyr::filter(symbol %in% top_genes) |> 
		ggplot(aes(y = symbol, x = .data[[x_var]], size = avg_log2FC, color = p_val_adj)) + 
		geom_point() +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

find_enrichment_by_cluster <- function(tbl_list){
	# browser()
	
	tbl_list |> 
		purrr::compact() |> 
		map(tibble::column_to_rownames, "symbol") |> 
		map(enrichment_analysis, annotate = FALSE) |>
		identity()
}



tar_load(c("cis_diffex_clones_for_each_cluster", "numbat_rds_files", "large_clone_comparisons", "cluster_orders", "cluster_dictionary", "debranched_seus_2p", "debranched_seus_6p", "debranched_seus_16q"))

# load seurat paths and objects ------------------------------

# 1q------------------------------
seus_1q <- dir_ls("output/seurat/integrated_1q/") %>%
	set_names(fs::path_file(.)) |> 
	# map(readRDS) |>
	identity()

integrated_seu_1q <- readRDS('output/seurat/integrated_1q/integrated_seu_1q_complete.rds')

# 2p ------------------------------
seus_2p <- debranched_seus_2p %>%
	set_names(fs::path_file(.)) |> 
	# map(readRDS) |>
	identity()

# 6p ------------------------------
seus_6p <- debranched_seus_6p |> 
	set_names(fs::path_file(.)) |> 
	# map(readRDS) |>
	identity()

# 16q ------------------------------
seus_16q <- dir_ls("output/seurat/integrated_16q/") %>%
	set_names(fs::path_file(.)) |> 
	# map(readRDS) |>
	identity()

integrated_seu_16q <- readRDS('output/seurat/integrated_16q/integrated_seu_16q_complete.rds')


# diffex per cluster ------------------------------

# 1q ------------------------------

all_diffex_by_cluster_1q <- map(seus_1q[1:4], find_diffex_bw_clones_for_each_cluster, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all")

saveRDS(all_diffex_by_cluster_1q, "results/1q_all_diffex_by_cluster_1q.rds")

all_diffex_by_cluster_1q <- readRDS("results/1q_all_diffex_by_cluster.rds")

all_diffex_by_cluster_1q <- 
	unlist(all_diffex_by_cluster_1q) |> 
	set_names(str_extract(unlist(all_diffex_by_cluster_1q), "SRR[0-9]*")) |> 
	map(read_csv)

all_diffex_by_cluster_1q <- map(all_diffex_by_cluster_1q, ~split(.x, .x$cluster)) |>
	purrr::list_transpose() |>
	map(~purrr::compact(.x)) |>
	identity()

cluster_diffex_plot_1q <- 
	all_diffex_by_cluster_1q |> 
	map(dplyr::bind_rows, .id = "sample_id") |> 
	map(dplyr::filter, p_val_adj <= 0.1) |>
	dplyr::bind_rows(.id = "clusters") |> 
	dplyr::group_by(symbol) %>%
	dplyr::filter(all_same_sign(avg_log2FC)) |> 
	dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
	dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
	dplyr::arrange(desc(abs_mean_FC)) |>
	dplyr::mutate(sample_cluster = glue("{sample_id}_{cluster}")) |> 
	select_recurrent_genes(x_var = "sample_cluster", 1) |>
	identity()

# 2p ------------------------------

all_diffex_by_cluster_2p <- map(seus_2p, find_diffex_bw_clones_for_each_cluster, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all")

saveRDS(all_diffex_by_cluster_2p, "results/all_diffex_by_cluster_2p.rds")

all_diffex_by_cluster_2p <- readRDS("results/all_diffex_by_cluster_2p.rds")

all_diffex_by_cluster_2p <- 
	unlist(all_diffex_by_cluster_2p) |> 
	set_names(str_extract(unlist(all_diffex_by_cluster_2p), "SRR[0-9]*")) |> 
	map(read_csv) |> 
	map(dplyr::bind_rows, .id = "clusters") |> 
	map(dplyr::filter, p_val_adj <= 0.1) |>
	dplyr::bind_rows(.id = "sample_id") |> 
	dplyr::group_by(symbol) %>%
	dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
	dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
	dplyr::arrange(desc(abs_mean_FC)) |>
	# purrr::discard(~ nrow(.x) == 0) |>
	identity()

cluster_diffex_plot_2p <- 
	all_diffex_by_cluster_2p |> 
	dplyr::mutate(sample_cluster = glue("{sample_id}_{cluster}")) |> 
	select_recurrent_genes(x_var = "sample_cluster", 2)

# 6p ------------------------------

all_diffex_by_cluster_6p <- map(seus_6p, find_diffex_bw_clones_for_each_cluster, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all")

saveRDS(all_diffex_by_cluster_6p, "results/all_diffex_by_cluster_6p.rds")

all_diffex_by_cluster_6p <- readRDS("results/all_diffex_by_cluster_6p")

all_diffex_by_cluster_6p <- 
	unlist(all_diffex_by_cluster_6p) |> 
	set_names(str_extract(unlist(all_diffex_by_cluster_6p), "SRR[0-9]*")) |> 
	map(read_csv)

# 16q ------------------------------





# debug(find_diffex_bw_clones_for_each_cluster)
# debug(clone_diff_per_cluster)
# debug(make_clone_comparison)

# cis_diffex <- map(seus[1:4], find_diffex_bw_clones_for_each_cluster, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "cis")

# diffex by clone ------------------------------

# 1q ------------------------------
all_diffex_by_clone_1q <- map(seus_1q[1:4], find_diffex_clones, numbat_rds_files, large_clone_comparisons, location = "all")

saveRDS(all_diffex_by_clone, "results/1q_all_diffex_by_clone.rds")

all_diffex_by_clone <- readRDS("results/1q_all_diffex_by_clone.rds")

pull_diffex_by_scna <- function(diffex_list, scna_of_interest = "1q\\+"){
	# browser()
	retained_diffexes <- str_detect(names(diffex_list), scna_of_interest)
	
	diffex_list <- diffex_list[retained_diffexes]
	
	return(diffex_list)
}

diffex_1q <- map(all_diffex_by_clone, pull_diffex_by_scna, scna_of_interest = "1q") |> 
	purrr::list_flatten() |> 
	dplyr::bind_rows(.id = "sample_id") |> 
	dplyr::filter(p_val_adj <= 0.05) |> 
	dplyr::group_by(symbol) %>%
	dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
	dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
	dplyr::arrange(desc(abs_mean_FC)) |>
	identity()

select_recurrent_genes(diffex_1q, x_var = "sample_id", 4) + 
	labs(title = "diffex by clone 1q")


# diffex by cluster ------------------------------

debug(find_diffex_bw_clones_for_each_cluster)

integrated_seu <- readRDS('output/seurat/integrated_1q/integrated_seu_1q_complete.rds')

cluster_seus <- SplitObject(integrated_seu, split.by = "clusters")

integrated_cluster_diffexes <- map(cluster_seus, FindMarkers, group.by = "scna", ident.1 = "w_scna", ident.2 = "wo_scna") |> 
	map(dplyr::filter, p_val_adj < 1) |> 
	map(tibble::rownames_to_column, "symbol") |> 
	# dplyr::bind_rows(.id = "clusters") |> 
	identity()

aggregated_expression <- AggregateExpression(integrated_seu, group.by = c("batch", "clusters", "scna"), return.seurat = TRUE)

aggregated_expression$group <- paste(aggregated_expression$clusters, aggregated_expression$scna, sep = "_")

# Idents(aggregated_expression) <- "group"

bulk.mono.de <- FindMarkers(object = aggregated_expression, 
														group.by = "scna",
														ident.1 = "w-scna",
														ident.2 = "wo-scna",
														test.use = "wilcox")
head(bulk.mono.de, n = 15)

FindMarkers(aggregated_expression, group.by = "")





enrichments_per_cluster <- 
	diffexes_by_cluster |> 
	map(find_enrichment_by_cluster)

combined_enrichment_plots <-
	enrichments_per_cluster |> 
	map(clusterProfiler::merge_result) |> 
	purrr::discard(~ nrow(.x) == 0) |> 
	map(plot_enrichment, result_slot = "compareClusterResult") |>
	imap(~{.x + labs(title = .y)}) |> 
	identity()


	


	
test1 <- 
	test0 |> 
	dplyr::bind_rows(.id = "cluster") |> 
	dplyr::group_by(symbol, cluster) %>%
	dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
	dplyr::mutate(recurrence = dplyr::n()) %>%
	dplyr::arrange(desc(recurrence), desc(abs_mean_FC)) |> 
	identity()

