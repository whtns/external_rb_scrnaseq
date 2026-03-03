# Plot Functions (103)

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @param group.by Parameter for group.by
#' @return Function result
#' @export
find_cluster_pairwise_distance <- function(seu, group.by) {
	df <- FetchData(seu, c(group.by, "S.Score", "G2M.Score")) |> 
		# dplyr::mutate(`G2M.Score` = (`G2M.Score` - min(`G2M.Score`)) / (max(`G2M.Score`) - min(`G2M.Score`))) |>
		# dplyr::mutate(`S.Score` = (`S.Score` - min(`S.Score`)) / (max(`S.Score`) - min(`S.Score`))) |>
		dplyr::rename(Cluster = .data[[group.by]]) |> 
		# dplyr::filter(Cluster %in% c("1", "3")) |> 
		dplyr::mutate(Cluster = as.character(Cluster))
	
	
	comparisons <- split(df, df$Cluster) |> 
		map(dplyr::select, -Cluster)
	
	cluster_pairs <- expand.grid(names(comparisons), names(comparisons)) |> 
		dplyr::filter(!Var1 == Var2) |> 
		dplyr::rename(x = Var1, y = Var2)
	
	
	# cluster_pairs <-
	# combn(names(comparisons), 2) |> 
	# 	t() |>
	# 	as_tibble() |>
	# 	dplyr::rename(x = V1, y = V2) |>
	# 	identity()
	
	avg_dists <- pmap_dbl(cluster_pairs, calculate_avg_cluster_distance, comparisons)
	
	cluster_pairs$avg_dist <- avg_dists
	
	cluster_order <- 
		cluster_pairs |> 
		dplyr::filter(avg_dist < 1) |> 
		dplyr::group_by(x) |>
		dplyr::summarize(mean_dist = mean(avg_dist)) |> 
		dplyr::arrange(mean_dist) |> 
		dplyr::pull(x) |> 
		as.character() |> 
		identity()
	
	# cluster_order <- unique(c(rbind(cluster_order$x, cluster_order$y)))
	# 
	cluster_pairs <- cluster_pairs |>
		dplyr::mutate(x = factor(x, levels = cluster_order))
	
	return(cluster_pairs)
}

#' Perform enrichment analysis
#'
#' @param enrichment_list Parameter for enrichment list
#' @param common_seus Parameter for common seus
#' @param plot_path File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
compare_corresponding_enrichments <- function(enrichment_list, common_seus = NULL, plot_path = tempfile(tmpdir = "results", fileext = ".pdf"), ...) {

	enrichment_results <- enrichment_list |> 
		map("enrichment") |> 
		map(1) |> 
		purrr::list_flatten() |>
		map(clusterProfiler::setReadable, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID") |>
		identity()
	
	enrichment_tables <-
		enrichment_results |>
		map(as_tibble)
	
	# tables_path <- writexl::write_xlsx(enrichment_tables, glue("results/corresponding_states_enrichment_{file_name}.xlsx"))
	
	enrichment_results <- 
		enrichment_results |> 
		map(dplyr::filter, p.adjust <= 0.1)
	
	common_terms <- enrichment_results[common_seus] |> 
		map(slot, "result") |> 
		map("ID") |>
		purrr::reduce(base::intersect) |> 
		identity()
	
	single_enrichment_plots <-  
		enrichment_results |> 
		map(plot_enrichment, showCategory = 20) |>
		imap(~ {
			
			# new_y_labels <- 
			# 	.x$data |> 
			# 	tibble::rownames_to_column("term") |> 
			# 	dplyr::mutate(label = str_remove(term, "^[a-z]*\\.")) |> 
			# 	dplyr::mutate(label = ifelse(label %in% common_terms, glue("**{label}**"), label)) |>
			# 	dplyr::pull(label) |> 
			# 	identity()
			
			.x + labs(title = .y) + 
				# scale_y_discrete(labels = new_y_labels) + 
				theme(
					axis.text.y = ggtext::element_markdown()
				)
		})
	
	merge_enrichment_plot <- enrichment_results |> 
		map(dplyr::filter, p.adjust <= 0.1) |> 
		clusterProfiler::merge_result() |> 
		plot_enrichment(p_val_cutoff = 0.05, result_slot = "compareClusterResult")

	
	pdf(plot_path, ...)
	print(single_enrichment_plots)
	print(merge_enrichment_plot)
	dev.off()
	
	return(plot_path)
}

#' Create a plot visualization
#'
#' @param sample_diffex_list Parameter for sample diffex list
#' @param filename File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_corresponding_enrichment <- function(sample_diffex_list, filename, ...){
	
	sample_id <- fs::path_file(filename)
	
	gene_sets <- map(c("H" = "H", "C6" = "C6"), ~msigdbr::msigdbr(species = "human", category = .x)) |> 
		map(dplyr::select, gs_name, entrez_gene) |> 
		dplyr::bind_rows()
	
	enrichments <- sample_diffex_list |> 
		map("all") |> 
		map(prep_for_enrichment, TERM2GENE = gene_sets)
	
	
#' Filter data based on specified criteria
#'
#' @param enrichment_result Parameter for enrichment result
#' @return Filtered data
#' @export
filter_enrichment_result <- function(enrichment_result){
		enrichment_result@result[enrichment_result@result$p.adjust < 0.25,] 
		return(enrichment_result)
	}
	
	enrichment_plots <- 
		enrichments |> 
		map(filter_enrichment_result) |> 
		map(plot_enrichment) |>
		imap(~{.x + labs(title = sample_id, subtitle = .y)}) |>
		map(~{.x + scale_y_discrete(labels = function(x) str_wrap(str_replace_all(str_remove(x, "HALLMARK"), "_", " "), width = 20))}) |>
		identity()
	
	enrichment_plots[[1]] + 
		scale_y_discrete(labels = function(x) str_wrap(str_replace_all(str_remove(x, "HALLMARK"), "_", " "), width = 10)) +
		# theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
		NULL
	
	enrichment_tables <-
		enrichments |> 
		map(clusterProfiler::setReadable, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID") |>
		map(as_tibble) |> 
		identity()
	
	tables_path <- writexl::write_xlsx(enrichment_tables, glue("results/corresponding_states_enrichment_{sample_id}.xlsx"))
	
	plot_path <- tempfile(tmpdir = "results", fileext = ".pdf")

	pdf(plot_path)
	print(enrichment_plots)
	dev.off()
	
	return(list("plot" = plot_path, "table" = tables_path, "enrichment" = enrichments))
	
}
#' Create a plot visualization
#'
#' @param sample_diffex_list Parameter for sample diffex list
#' @param filename File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_fig_03_09 <- function(sample_diffex_list, filename, ...){
	
	sample_id <- fs::path_file(filename)
	
	gene_sets <- map(c("H" = "H", "C6" = "C6"), ~msigdbr::msigdbr(species = "human", category = .x)) |> 
		map(dplyr::select, gs_name, entrez_gene) |> 
		dplyr::bind_rows()
	
	enrichments <- sample_diffex_list |> 
		map("all") |> 
		map(prep_for_enrichment, TERM2GENE = gene_sets)
	
	enrichment_plots <- 
		enrichments |> 
		map(filter_enrichment_result) |> 
		map(plot_enrichment) |>
		imap(~{.x + labs(title = sample_id, subtitle = .y)}) |>
		identity()
	
	enrichment_tables <-
		enrichments |> 
		map(clusterProfiler::setReadable, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID") |>
		map(as_tibble) |> 
		identity()
	
	tables_path <- writexl::write_xlsx(enrichment_tables, glue("results/corresponding_states_enrichment_{sample_id}.xlsx"))
	
	plot_path <- glue("results/{sample_id}_corresponding_clusters_diffex_enrichment.pdf", h = 10)
	pdf(plot_path)
	print(enrichment_plots)
	dev.off()
	
	return(list("plot" = plot_path, "table" = tables_path, "enrichment" = enrichments))
	
}

