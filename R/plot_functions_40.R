# Plot Functions (139)
#' Create a plot visualization
#'
#' @param unfiltered_seu_path File path
#' @param filtered_seu_path File path
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

plot_effect_of_filtering_old <- function(unfiltered_seu_path, filtered_seu_path, group.by = "SCT_snn_res.0.6") {
  
	#
	
	plot_list <- list()
	
	sample_id <- str_extract(unfiltered_seu_path, "SRR[0-9]*")
	
	unfiltered_seu <- unfiltered_seu_path
	
	filtered_seu <- filtered_seu_path
	
	fs::dir_create("results/effect_of_filtering")
	
	# distribution ------------------------------
	dir_create("results/effect_of_filtering/distribution")
	plot_distribution_of_clones_across_clusters(filtered_seu, sample_id, var_x = "scna", var_y = group.by)
	fs::dir_create(glue("results/effect_of_filtering/distribution/filtered/"))
	plot_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_filtered_distribution.pdf")
	plot_list["filtered_distribution"] <- plot_path
	ggsave(plot_path, height = 4, width = 8)
	
	filtered_dist_tables <- table_distribution_of_clones_across_clusters(filtered_seu, sample_id, clusters = group.by)
	
	table_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_filtered_distribution.xlsx")
	plot_list["filtered_distribution_tables"] <- table_path
	writexl::write_xlsx(filtered_dist_tables, table_path)
	
	plot_distribution_of_clones_across_clusters(unfiltered_seu, sample_id, var_y = group.by)
	fs::dir_create(glue("results/effect_of_filtering/distribution/filtered"))
	plot_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_unfiltered_distribution.pdf")
	plot_list["unfiltered_distribution"] <- plot_path
	ggsave(plot_path, height = 4, width = 8)
	
	unfiltered_dist_tables <- table_distribution_of_clones_across_clusters(unfiltered_seu, sample_id, clusters = group.by)
	
	table_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_unfiltered_distribution.xlsx")
	plot_list["unfiltered_distribution_tables"] <- table_path
	writexl::write_xlsx(unfiltered_dist_tables, table_path)
	
	# abbreviation markers ------------------------------
	dir_create("results/effect_of_filtering/abbreviation")
	(plot_markers(filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
			labs(title = "filtered")) +
		(plot_markers(unfiltered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
		 	labs(title = "unfiltered")) +
		plot_annotation(title = sample_id)
	
	plot_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_markers.pdf")
	plot_list["abbreviation_markers"] <- plot_path
	ggsave(plot_path, height = 12, width = 15)
	
	filtered_marker_tables <- table_cluster_markers(filtered_seu)
	
	table_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_filtered_markers.xlsx")
	plot_list["filtered_marker_tables"] <- table_path
	writexl::write_xlsx(filtered_marker_tables, table_path)
	
	unfiltered_marker_tables <- table_cluster_markers(unfiltered_seu) %>%
		purrr::compact()
	
	table_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_unfiltered_markers.xlsx")
	plot_list["unfiltered_marker_tables"] <- table_path
	writexl::write_xlsx(unfiltered_marker_tables, table_path)
	
	heatmap_features <-
		table_cluster_markers(unfiltered_seu) %>%
		pluck(group.by) %>%
		group_by(Cluster) %>%
		slice_head(n = 10) %>%
		dplyr::pull(Gene.Name) %>%
		identity()
	
	ggplotify::as.ggplot(
		seu_complex_heatmap(unfiltered_seu,
												features = heatmap_features,
												group.by = c(group.by, "Phase", "scna"),
												col_arrangement = c(group.by, "Phase", "scna"),
												cluster_rows = FALSE, use_raster = TRUE
		)
	) +
		labs(title = sample_id)
	
	plot_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_heatmap.pdf")
	plot_list["abbreviation_heatmap"] <- plot_path
	ggsave(plot_path, height = 8, width = 8)
	
	# nCount_gene umaps ------------------------------
	
	unfiltered_seu$log_nCount_gene <- log1p(unfiltered_seu$nCount_gene)
	filtered_seu$log_nCount_gene <- log1p(filtered_seu$nCount_gene)
	
	(FeaturePlot(unfiltered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
			labs(title = "unfiltered")) +
		(FeaturePlot(filtered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
		 	labs(title = "filtered")) +
		plot_annotation(title = sample_id)
	
	fs::dir_create(glue("results/effect_of_filtering/nCount_gene"))
	plot_path <- glue("results/effect_of_filtering/nCount_gene/{sample_id}_nCount_gene_umaps.pdf")
	plot_list["nCount_gene_umaps"] <- plot_path
	ggsave(plot_path, height = 8, width = 10)
	
	# abbreviation umaps ------------------------------
	mycols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["abbreviation"]])))
	
	(DimPlot(unfiltered_seu, group.by = "abbreviation", cols = mycols) +
			labs(title = "unfiltered")) +
		(DimPlot(filtered_seu, group.by = "abbreviation", cols = mycols) +
		 	labs(title = "filtered")) +
		# (DimPlot(regressed_seu, group.by = "gene_snn_res.0.2") +
		#   labs(title = "regressed")) +
		plot_annotation(title = sample_id)
	
	fs::dir_create(glue("results/effect_of_filtering/abbreviation"))
	plot_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_umaps.pdf")
	plot_list["abbreviation_umaps"] <- plot_path
	ggsave(plot_path, height = 8, width = 10)
	
	# scna umaps ------------------------------
	mycols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["scna"]])))
	
	unfiltered_seu@meta.data$scna <- vec_split_label_line(unfiltered_seu@meta.data$scna, 3)
	filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)
	
	(DimPlot(unfiltered_seu, group.by = "scna", cols = mycols) +
			labs(title = "unfiltered")) +
		(DimPlot(filtered_seu, group.by = "scna", cols = mycols) +
		 	labs(title = "filtered")) +
		# (DimPlot(regressed_seu, group.by = "gene_snn_res.0.2") +
		#     labs(title = "regressed")) +
		plot_annotation(title = sample_id)
	
	fs::dir_create(glue("results/effect_of_filtering/scna"))
	plot_path <- glue("results/effect_of_filtering/scna/{sample_id}_scna_umaps.pdf")
	plot_list["scna_umaps"] <- plot_path
	ggsave(plot_path, height = 8, width = 10)
	
	return(plot_list)
}

#' Create a plot visualization
#'
#' @param unfiltered_seu_path File path
#' @param filtered_seu_path File path
#' @param group.by Character string (default: "gene_snn_res.0.2")
#' @param cluster_dictionary Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_effect_of_filtering <- function(unfiltered_seu_path, filtered_seu_path = NULL, group.by = "gene_snn_res.0.2", cluster_dictionary, plot_path = "results/fig_02.01.1.pdf") {
  
  
  
  
  
  
	
	sample_id <- str_extract(unfiltered_seu_path, "SRR[0-9]*")
	
	unfiltered_seu <- unfiltered_seu_path
	
	original_filtered_seu <- filtered_seu_path

  plot_list <- list()
  
  removed_clusters <- cluster_dictionary[[sample_id]] |> 
  	dplyr::filter(remove == 1) |> 
  	dplyr::pull(`gene_snn_res.0.2`)

  # removed_clusters <- removed_clusters[removed_clusters$sample_id == sample_id,][["gene_snn_res.0.2"]]
  
  filtered_seu <- unfiltered_seu[,!unfiltered_seu$gene_snn_res.0.2 %in% removed_clusters]

  # abbreviation markers ------------------------------
  (plot_markers(unfiltered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 5) +
  		labs(title = "unfiltered") + 
   	theme(legend.position = "none") + 
   	NULL) +
  	(plot_markers(original_filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 5) +
    labs(title = "filtered") + 
    	theme(axis.title.y = element_blank()) +
    	NULL) +
    plot_annotation(title = sample_id) + 
  	plot_layout(guides='collect')

  tmp_path <- tempfile(fileext = ".pdf")
  plot_list["abbreviation_markers"] <- tmp_path
  ggsave(tmp_path, height = 6, width = 9)

  # scna umaps ------------------------------
  scna_cols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["scna"]])))
  
  unfiltered_seu@meta.data$scna <- vec_split_label_line(unfiltered_seu@meta.data$scna, 3)
  filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)

  # abbreviation umaps ------------------------------
  col_numbers <- sort(as.numeric(unique(unfiltered_seu@meta.data[[group.by]])))
  group_cols <- scales::hue_pal()(length(col_numbers)) |> 
  	set_names(col_numbers)

  tmp_path <- tempfile(fileext = ".pdf")
  plot_list["abbreviation_umaps"] <- tmp_path
  pdf(tmp_path, height = 3, width = 4.5)
  
  print(DimPlot(unfiltered_seu, group.by = "scna", cols = scna_cols) +
  		labs(title = "unfiltered") +
  		theme(axis.title.x = element_blank()))
  	
  	print(DimPlot(filtered_seu, group.by = "scna", cols = scna_cols) +
  	 	labs(title = "filtered") +
  	 	theme(axis.title.x = element_blank(),
  	 				axis.title.y = element_blank()))
  	
  	print(DimPlot(original_filtered_seu, group.by = "scna", cols = scna_cols) +
  	 	labs(title = "re-clustered") +
  	 	theme(axis.title.x = element_blank(),
  	 				axis.title.y = element_blank()))

  	print(DimPlot(unfiltered_seu, group.by = group.by, cols = group_cols))
   
  	print(DimPlot(filtered_seu, group.by = group.by, cols = group_cols) +
     	theme(axis.title.y = element_blank(),
     				legend.position = "none"))

  	print(DimPlot(original_filtered_seu, group.by = group.by, cols = group_cols) +
  	 	theme(axis.title.y = element_blank(),
  	 				legend.position = "none"))
  	
  	dev.off()
  
  plot_path <- qpdf::pdf_combine(plot_list, plot_path)

  return(plot_path)
}

#' Create a plot visualization
#'
#' @param filtered_seu_path File path
#' @param regressed_seu_path File path
#' @param resolution Parameter for resolution
#' @param filter_dropped_cluster Cluster information
#' @param regress_dropped_cluster Cluster information
#' @param n_features Parameter for n features
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_effect_of_regression <- function(filtered_seu_path, regressed_seu_path, resolution = 0.4, filter_dropped_cluster = NULL, regress_dropped_cluster = NULL, n_features = 2, ...) {
  
  
  
  
  
  
  #

  sample_id <- str_extract(filtered_seu_path, "SRR[0-9]*")

  regressed_seu <- readRDS(regressed_seu_path)


  regressed_seu$scna[regressed_seu$scna == ""] <- ".diploid"
  regressed_seu$scna <- factor(regressed_seu$scna)

  regressed_meta <- regressed_seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, starts_with("SCT_snn_res")) %>%
    tibble::column_to_rownames("cell") %>%
    dplyr::rename_with(
      ~ str_replace(.x, "SCT_snn_res.", "regressed."),
      starts_with("SCT_")
    ) %>%
    identity()

  filtered_seu <- readRDS(filtered_seu_path)

  filtered_seu$scna[filtered_seu$scna == ""] <- ".diploid"
  filtered_seu$scna <- factor(filtered_seu$scna)

  filtered_seu <- AddMetaData(filtered_seu, regressed_meta)

  filtered_seu <- find_all_markers(filtered_seu, colnames(regressed_meta))

  regressed_cluster <- glue("regressed.{resolution}")
  
  # if(!is.null(regress_dropped_cluster)){
  # 	regressed_seu <- regressed_seu[,!regressed_seu@meta.data[[regressed_cluster]] %in% regress_dropped_cluster]
  # }

  filtered_cluster <- glue("SCT_snn_res.{resolution}")
  
  # browser()
  if(!is.null(filter_dropped_cluster)){
  	filtered_seu <- filtered_seu[,!filtered_seu@meta.data[[filtered_cluster]] %in% filter_dropped_cluster]
  }

  regressed_features <-
    filtered_seu@misc$markers[[regressed_cluster]]$presto %>%
    group_by(Cluster) %>%
    slice_head(n = n_features) %>%
    dplyr::select(regressed_cluster = Cluster, Gene.Name) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
    dplyr::filter(Gene.Name %in% rownames(GetAssayData(filtered_seu, "SCT", "scale.data"))) %>%
    identity()

  filtered_features <-
    filtered_seu@misc$markers[[filtered_cluster]]$presto %>%
    group_by(Cluster) %>%
    slice_head(n = n_features) %>%
    dplyr::select(Cluster, Gene.Name) %>%
    dplyr::mutate(filtered_cluster = Cluster) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
    dplyr::filter(Gene.Name %in% rownames(GetAssayData(filtered_seu, "SCT", "scale.data"))) %>%
    identity()

  # filtered ------------------------------

  heatmap_features <- filtered_features %>%
    dplyr::left_join(regressed_features, by = "Gene.Name") %>%
    dplyr::mutate(regressed_cluster = replace_na(regressed_cluster, "NA")) %>%
    identity()

  filtered_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(filtered_seu,
      features = heatmap_features$Gene.Name,
      group.by = c("G2M.Score", "S.Score", "scna", regressed_cluster),
      col_arrangement = c(filtered_cluster, "scna"),
      cluster_rows = FALSE,
      column_split = sort(filtered_seu@meta.data[[filtered_cluster]]),
      row_split = rev(heatmap_features$filtered_cluster),
      row_title_rot = 0,
      use_raster = TRUE,
      # right_annotation = row_ha
    )
  ) +
    labs(title = "filtered") +
    theme()

  cc_data <- FetchData(filtered_seu, c(filtered_cluster, regressed_cluster, "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(.data[[filtered_cluster]]) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    # dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate({{ filtered_cluster }} := as.factor(.data[[filtered_cluster]])) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  filtered_facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[filtered_cluster]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[filtered_cluster]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    facet_wrap(~ .data[[filtered_cluster]], ncol = 3) +
    theme_light() +
    # geom_label(data = labels,
    # 					 aes(label = label),
    # 					 # x = Inf,
    # 					 # y = -Inf,
    # 					 x = max(cc_data$S.Score)+0.05,
    # 					 y = max(cc_data$G2M.Score)-0.1,
    # 					 hjust=1,
    # 					 vjust=1,
    # 					 inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(title = "filtered") +
    # guides(color = "none") +
    NULL

  # regressed ------------------------------

  heatmap_features <- regressed_features %>%
    dplyr::left_join(filtered_features, by = "Gene.Name") %>%
    dplyr::mutate(filtered_cluster = replace_na(filtered_cluster, "NA")) %>%
    identity()

  regressed_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(filtered_seu,
      features = heatmap_features$Gene.Name,
      group.by = c("G2M.Score", "S.Score", "scna", filtered_cluster),
      col_arrangement = c(regressed_cluster, "scna"),
      cluster_rows = FALSE,
      column_split = sort(filtered_seu@meta.data[[regressed_cluster]]),
      row_split = rev(heatmap_features$regressed_cluster),
      row_title_rot = 0,
      use_raster = TRUE
      # right_annotation = row_ha
    )
  ) +
    labs(title = "regressed") +
    theme()

  cc_data <- FetchData(filtered_seu, c(filtered_cluster, regressed_cluster, "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(.data[[regressed_cluster]]) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    # dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate({{ regressed_cluster }} := as.factor(.data[[regressed_cluster]])) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  regressed_facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[regressed_cluster]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[regressed_cluster]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    facet_wrap(~ .data[[regressed_cluster]], ncol = 3) +
    theme_light() +
    # geom_label(data = labels,
    # 					 aes(label = label),
    # 					 # x = Inf,
    # 					 # y = -Inf,
    # 					 x = max(cc_data$S.Score)+0.05,
    # 					 y = max(cc_data$G2M.Score)-0.1,
    # 					 hjust=1,
    # 					 vjust=1,
    # 					 inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(title = "regressed") +
    # guides(color = "none") +
    NULL


  # patchworks ------------------------------
  
  plot_list = list(
  	"A" = filtered_heatmap,
  	"B" = regressed_heatmap,
  	"C" = filtered_facet_cell_cycle_plot,
  	"D" = regressed_facet_cell_cycle_plot,
  	"E" = plot_spacer(),
  	"F" = plot_spacer()
  )
  
  layout <- "
              AAAABBBB
              AAAABBBB
              AAAABBBB
              CCCEDDDF
              CCCEDDDF
              "
  
  mypatch <- wrap_plots(
  	plot_list
  ) +
  	plot_layout(design = layout) +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/{sample_id}_regression_effects.pdf")
  ggsave(plot_path, plot = mypatch, ...)

  # markerplot <- plot_seu_marker_heatmap(filtered_seu_path, nb_path = numbat_rds_file, clone_simplifications = large_clone_simplifications, ...)
  return(plot_path)
}

#' Create a plot visualization
#'
#' @param filtered_seu_path File path
#' @param regressed_seu_path File path
#' @param resolution Character string (default: "0.4")
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return ggplot2 plot object
#' @export
plot_effect_of_regression_old <- function(filtered_seu_path, regressed_seu_path, resolution = "0.4", group.by = "SCT_snn_res.0.6") {
  
  
  
  
  
  
  #

  plot_list <- list()

  sample_id <- str_extract(filtered_seu_path, "SRR[0-9]*")

  filtered_seu <- filtered_seu_path

  regressed_seu <- regressed_seu_path

  fs::dir_create("results/effect_of_regression")

  # percent.mt umaps ------------------------------

  (FeaturePlot(filtered_seu, features = "percent.mt") +
    labs(title = "filtered")) +
    (FeaturePlot(regressed_seu, features = "percent.mt") +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/percent_mt"))
  plot_path <- glue("results/effect_of_regression/percent_mt/{sample_id}_percent_mt_umaps.pdf")
  plot_list["percent_mt_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)

  # distribution ------------------------------

  plot_distribution_of_clones_across_clusters(filtered_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_regression/distribution/filtered/"))
  plot_path <- glue("results/effect_of_regression/distribution/filtered/{sample_id}_filtered_distribution.pdf")
  plot_list["filtered_distribution"] <- plot_path
  ggsave(plot_path, height = 4, width = 8)

  filtered_dist_tables <- table_distribution_of_clones_across_clusters(filtered_seu, sample_id, clusters = group.by)

  table_path <- glue("results/effect_of_regression/distribution/filtered/{sample_id}_filtered_distribution.xlsx")
  plot_list["filtered_distribution_tables"] <- table_path
  writexl::write_xlsx(filtered_dist_tables, table_path)

  plot_distribution_of_clones_across_clusters(regressed_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_regression/distribution/regression"))
  plot_path <- glue("results/effect_of_regression/distribution/regression/{sample_id}_regressed_distribution.pdf")
  plot_list["regressed_distribution"] <- plot_path
  ggsave(plot_path, height = 4, width = 8)

  regressed_dist_tables <- table_distribution_of_clones_across_clusters(regressed_seu, sample_id, clusters = group.by)

  table_path <- glue("results/effect_of_regression/distribution/regression/{sample_id}_regressed_distribution.xlsx")
  plot_list["regressed_distribution_tables"] <- table_path
  writexl::write_xlsx(regressed_dist_tables, table_path)

  # abbreviation markers ------------------------------
  (plot_markers(filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
    labs(title = "filtered")) +
    (plot_markers(regressed_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_markers.pdf")
  plot_list["abbreviation_markers"] <- plot_path
  ggsave(plot_path, height = 12, width = 15)

  filtered_marker_tables <- table_cluster_markers(filtered_seu)

  table_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_filtered_markers.xlsx")
  plot_list["filtered_marker_tables"] <- table_path
  writexl::write_xlsx(filtered_marker_tables, table_path)

  regressed_marker_tables <- table_cluster_markers(regressed_seu) %>%
    purrr::compact()

  table_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_regressed_markers.xlsx")
  plot_list["regressed_marker_tables"] <- table_path
  writexl::write_xlsx(regressed_marker_tables, table_path)

  heatmap_features <-
    table_cluster_markers(regressed_seu) %>%
    pluck(group.by) %>%
    group_by(Cluster) %>%
    slice_head(n = 10) %>%
    dplyr::pull(Gene.Name) %>%
    identity()

  ggplotify::as.ggplot(
    seu_complex_heatmap(regressed_seu,
      features = heatmap_features,
      group.by = c(group.by, "Phase", "scna"),
      col_arrangement = c(group.by, "Phase", "scna"),
      cluster_rows = FALSE
    )
  ) +
    labs(title = sample_id)

  plot_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_heatmap.pdf")
  plot_list["abbreviation_heatmap"] <- plot_path
  ggsave(plot_path, height = 8, width = 8)

  # # nCount_gene umaps ------------------------------
  #
  # filtered_seu$log_nCount_gene <- log1p(filtered_seu$nCount_gene)
  # regressed_seu$log_nCount_gene <- log1p(regressed_seu$nCount_gene)
  #
  # (FeaturePlot(filtered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
  #     labs(title = "filtered")) +
  #   (FeaturePlot(regressed_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
  #      labs(title = "regressed")) +
  #   plot_annotation(title = sample_id)
  #
  # fs::dir_create(glue("results/effect_of_regression/nCount_gene"))
  # plot_path = glue("results/effect_of_regression/nCount_gene/{sample_id}_nCount_gene_umaps.pdf")
  # plot_list["nCount_gene_umaps"] = plot_path
  # ggsave(plot_path, height = 8, width = 10)

  # # cluster umaps ------------------------------
  # mycols = scales::hue_pal()(length(unique(filtered_seu@meta.data[[group.by]])))
  #
  # (DimPlot(filtered_seu, group.by = group.by, cols = mycols) +
  #   labs(title = "filtered")) +
  # (DimPlot(regressed_seu, group.by = group.by, cols = mycols) +
  #   labs(title = "regressed")) +
  #   plot_annotation(title = sample_id)
  #
  # fs::dir_create(glue("results/effect_of_regression/cluster"))
  # plot_path = glue("results/effect_of_regression/cluster/{sample_id}_cluster_umaps.pdf")
  # plot_list["cluster_umaps"] = plot_path
  # ggsave(plot_path, height = 8, width = 10)
  #
  # abbreviation umaps ------------------------------
  mycols <- scales::hue_pal()(length(unique(filtered_seu@meta.data[[group.by]])) + 3)

  (DimPlot(filtered_seu, group.by = group.by, cols = mycols) +
    labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = group.by, cols = mycols) +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/abbreviation"))
  plot_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_umaps.pdf")
  plot_list["abbreviation_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)

  # scna umaps ------------------------------
  mycols <- scales::hue_pal()(length(unique(filtered_seu@meta.data[["scna"]])) + 3)

  # filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)
  # regressed_seu@meta.data$scna <- vec_split_label_line(regressed_seu@meta.data$scna, 3)

  (DimPlot(filtered_seu, group.by = "scna", cols = mycols) +
    labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = "scna", cols = mycols) +
      labs(title = "regressed")) +
    # (DimPlot(regressed_filtered_seu, group.by = group.by) +
    #     labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/scna"))
  plot_path <- glue("results/effect_of_regression/scna/{sample_id}_scna_umaps.pdf")
  plot_list["scna_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)

  # # scna markers ------------------------------
  # regressed_seu$scna[regressed_seu$scna == ""] <- "none"
  #
  # plot_markers(regressed_seu, "scna", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10, unique_markers = TRUE) +
  #   labs(title = sample_id)
  #
  # plot_path = glue("results/effect_of_regression/{sample_id}_scna_markers.pdf")
  # plot_list["scna_markers"] = plot_path
  # ggsave(plot_path, height = 8, width = 6)

  # browseURL(glue("results/effect_of_regression/{sample_id}_markers_by_scna.pdf"))
  #
  # #   original_mt_plot <-
  # #     FeaturePlot(regressed_seu, features = "percent.mt") +
  # #     labs(title = "filtered")
  # #
  # #   regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "percent.mt") +
  # #   labs(title = "regressed")
  # #
  # #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  # #     plot_annotation(title = sample_id)
  # #
  # # ggsave(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"), heigh = 6, width = 10)
  # # browseURL(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"))
  #
  # regressed_seu <- AddModuleScore(regressed_seu, subtype_markers)
  # regressed_seu0 <- AddModuleScore(regressed_seu0, subtype_markers)
  #
  #   original_mt_plot <-
  #     FeaturePlot(regressed_seu, features = "Cluster1") +
  #     labs(title = "filtered")
  #
  #   regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "Cluster1") +
  #   labs(title = "regressed")
  #
  #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #     plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"))
  #
  # original_mt_plot <-
  #   FeaturePlot(regressed_seu, features = "Cluster2") +
  #   labs(title = "filtered")
  #
  # regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "Cluster2") +
  #   labs(title = "regressed")
  #
  # wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"))

  # DimPlot(regressed_seu, group.by = group.by) +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = group.by) +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_louvain.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_louvain.pdf"))
  #
  # DimPlot(regressed_seu, group.by = "scna") +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = "scna") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_scna.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_scna.pdf"))
  #
  # DimPlot(regressed_seu, group.by = "Phase") +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = "Phase") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_phase.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_phase.pdf"))

  # phase umaps ------------------------------

  (DimPlot(filtered_seu, group.by = "Phase") +
    labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = "Phase") +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/phase"))
  plot_path <- glue("results/effect_of_regression/phase/{sample_id}_phase_umaps.pdf")
  plot_list["phase_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)


  return(plot_list)
}

