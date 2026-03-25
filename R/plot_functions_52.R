# Plot Functions (151)

#' Perform print table tally operation
#'
#' @param comparison_id Parameter for comparison id
#' @param myvec Parameter for myvec
#' @return Function result
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

print_table_tally <- function(comparison_id, myvec) {
  table_tally <- table(myvec)
  message(comparison_id)
  message(glue("{names(table_tally)}: {table_tally}; "))
}

#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param metavar Character string (default: "batch")
#' @param num_markers Parameter for num markers
#' @param selected_values Parameter for selected values
#' @param return_plotly Logical flag (default: FALSE)
#' @param marker_method Character string (default: "presto")
#' @param seurat_assay Character string (default: "gene")
#' @param hide_technical Parameter for hide technical
#' @param unique_markers Logical flag (default: FALSE)
#' @param p_val_cutoff Threshold value for filtering
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
poster_plot_markers <- function(seu, metavar = "batch", num_markers = 5, selected_values = NULL,
                                return_plotly = FALSE, marker_method = "presto", seurat_assay = "gene",
                                hide_technical = NULL, unique_markers = FALSE, p_val_cutoff = 1,
                                ...) {
  Idents(seu) <- seu[[]][[metavar]]
  seu <- find_all_markers(seu, metavar,
    seurat_assay = seurat_assay,
    p_val_cutoff = p_val_cutoff
  )
  marker_table <- seu@misc$markers[[metavar]][[marker_method]]
  markers <- marker_table %>%
    seuratTools::enframe_markers() %>%
    dplyr::mutate(dplyr::across(.fns = as.character))
  if (!is.null(hide_technical)) {
    markers <- purrr::map(markers, c)
    if (hide_technical == "pseudo") {
      markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
    } else if (hide_technical == "mito_ribo") {
      markers <- purrr::map(markers, ~ .x[!str_detect(
        .x,
        "^MT-"
      )])
      markers <- purrr::map(markers, ~ .x[!str_detect(
        .x,
        "^RPS"
      )])
      markers <- purrr::map(markers, ~ .x[!str_detect(
        .x,
        "^RPL"
      )])
    } else if (hide_technical == "all") {
      markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
      markers <- purrr::map(markers, ~ .x[!str_detect(
        .x,
        "^MT-"
      )])
      markers <- purrr::map(markers, ~ .x[!str_detect(
        .x,
        "^RPS"
      )])
      markers <- purrr::map(markers, ~ .x[!str_detect(
        .x,
        "^RPL"
      )])
    }
    min_length <- min(purrr::map_int(markers, length))
    markers <- purrr::map(markers, head, min_length) %>%
      dplyr::bind_cols()
  }
  if (unique_markers) {
    markers <- markers %>%
      dplyr::mutate(precedence = row_number()) %>%
      pivot_longer(-precedence, names_to = "group", values_to = "markers") %>%
      dplyr::arrange(markers, precedence) %>%
      dplyr::group_by(markers) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::arrange(
        group,
        precedence
      ) %>%
      tidyr::drop_na() %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(precedence = row_number()) %>%
      tidyr::pivot_wider(
        names_from = "group",
        values_from = "markers"
      ) %>%
      dplyr::select(-precedence)
  }
  sliced_markers <- markers %>%
    dplyr::slice_head(n = num_markers) %>%
    tidyr::pivot_longer(everything(),
      names_to = "group",
      values_to = "feature"
    ) %>%
    dplyr::arrange(group) %>%
    dplyr::distinct(feature, .keep_all = TRUE) %>%
    identity()
  if (!is.null(selected_values)) {
    seu <- seu[, Idents(seu) %in% selected_values]
    sliced_markers <- sliced_markers %>%
      dplyr::filter(group %in%
        selected_values) %>%
      dplyr::distinct(feature, .keep_all = TRUE)
  }
  vline_coords <- head(cumsum(table(sliced_markers$group)) +
    0.5, -1)
  sliced_markers <- dplyr::pull(sliced_markers, feature)
  seu[[metavar]][is.na(seu[[metavar]])] <- "NA"
  Idents(seu) <- metavar
  markerplot <- DotPlot(seu,
    assay = seurat_assay, features = sliced_markers,
    group.by = metavar, dot.scale = 3
  ) + ggplot2::theme(axis.text.x = ggplot2::element_text(
    size = 10,
    angle = 45, vjust = 1, hjust = 1
  ), axis.text.y = ggplot2::element_text(size = 10)) +
    ggplot2::scale_y_discrete(position = "left") + ggplot2::scale_x_discrete(limits = sliced_markers) +
    ggplot2::geom_vline(xintercept = vline_coords, linetype = 2) +
    NULL
  if (return_plotly == FALSE) {
    return(markerplot)
  }
  plot_height <- (150 * num_markers)
  plot_width <- (100 * length(levels(Idents(seu))))
  markerplot <- plotly::ggplotly(markerplot,
    height = plot_height,
    width = plot_width
  ) %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    identity()
  return(list(plot = markerplot, markers = marker_table))
}

#' Perform calc silhouette operation
#'
#' @param seu_path File path
#' @param assay Character string (default: "SCT")
#' @param reduction Character string (default: "pca")
#' @param dims Parameter for dims
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
calc_silhouette <- function(seu_path, assay = "SCT", reduction = "pca", dims = 1:30, ...) {
  #
  tumor_id <- str_extract(seu_path, "SRR[0-9]*")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  plot_path <- glue("results/{sample_id}_silhouette.pdf")

  seu <- readRDS(seu_path)

  resolutions <- glue("{assay}_snn_res.{seq(0.2, 1.0, by = 0.2)}") %>%
    set_names(.)

  for (resolution in resolutions) {
    dist.matrix <- dist(x = Embeddings(object = seu[[reduction]])[, dims])
    clusters <- seu@meta.data[[resolution]]
    sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)

    sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")

    seu[[sil_res]] <- sil[, 3]
  }
  plot_list <- vector(mode = "list", length = length(resolutions)) |>
    set_names(resolutions)

  metrics <- vector(mode = "list", length = length(resolutions)) |>
    set_names(resolutions)

  for (resolution in resolutions) {
    #

    mod_out <- utils::capture.output(Seurat::FindClusters(seu, resolution = as.numeric(str_replace(resolution, "SCT_snn_res.", ""))), type = "output")
    modularity <- as.numeric(unlist(strsplit(x = mod_out[7], split = ":"))[2])

    sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")

    seu_meta <-
      seu@meta.data |>
      tibble::rownames_to_column("cell") |>
      dplyr::arrange(.data[[resolution]], desc(.data[[sil_res]])) |>
      dplyr::mutate(cell = factor(cell, levels = unique(cell)))

    mean_sil <- mean(seu_meta[[sil_res]])

    plot_list[[resolution]] <- ggplot(data = seu_meta, aes(x = .data[["cell"]], y = .data[[sil_res]], fill = .data[[resolution]])) +
      geom_col(outlier.size = 0.1) +
      geom_hline(aes(yintercept = mean_sil)) +
      xlab("Method") +
      ylab("Silhoutte Metric") +
      labs(title = glue("{sample_id} {resolution}"), subtitle = glue("silhouette mean: {mean_sil}; modularity: {modularity}")) +
      NULL

    metrics[[resolution]] <- c(
      "resolution" = as.numeric(str_replace(resolution, "SCT_snn_res.", "")),
      "silhouette_mean" = mean_sil
    )
  }

  metrics_plot <- dplyr::bind_rows(metrics) |>
    ggplot(aes(x = resolution, y = silhouette_mean)) +
    geom_point() +
    geom_line() +
    labs(title = sample_id)

  pdf(plot_path, ...)
  print(metrics_plot)
  print(plot_list)
  dev.off()

  return(plot_path)
}

#' Perform find genes by arm operation
#'
#' @param gene_list Gene names or identifiers
#' @return Function result
#' @export
find_genes_by_arm <- function(gene_list) {
	cc_genes <- gene_list %>%
		tibble::enframe("phase_of_gene", "symbol") %>%
		tidyr::unnest(symbol) |>
		dplyr::left_join(annotables::grch38, by = "symbol") |>
		dplyr::mutate(seqnames = chr) |>
		dplyr::filter(!is.na(start)) |>
		as_granges() |>
		identity()
	
	get_arms_ranges() |>
		join_overlap_intersect(cc_genes) |>
		as_tibble() |>
		dplyr::mutate(seqnames = str_pad(seqnames, 2, pad = "0")) |>
		dplyr::arrange(seqnames, arm) |>
		dplyr::select(seqnames, arm, everything())
}

#' Perform find cc genes by arm operation
#'
#' @return Function result
#' @export
find_cc_genes_by_arm <- function() {
  cc_genes <- Seurat::cc.genes.updated.2019 %>%
    tibble::enframe("phase_of_gene", "symbol") %>%
    tidyr::unnest(symbol) |>
    dplyr::left_join(annotables::grch38, by = "symbol") |>
    dplyr::mutate(seqnames = chr) |>
    dplyr::filter(!is.na(start)) |>
    as_granges() |>
    identity()

  get_arms_ranges() |>
    join_overlap_intersect(cc_genes) |>
    as_tibble() |>
    dplyr::mutate(seqnames = str_pad(seqnames, 2, pad = "0")) |>
    dplyr::arrange(seqnames, arm) |>
    dplyr::select(seqnames, arm, everything())
}

plot_feature_in_seu <- function(numbat_rds_file, ...) {
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[, !is.na(seu$clone_opt)]

  Seurat::FeaturePlot(seu, ...)
}

plot_clone_pearls <- function(seu_path, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), var_y = "clusters") {
	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
  cluster_order <- read_cluster_orders_table(file_id = file_id)
  if (length(cluster_order) > 0) {
    cluster_order <- cluster_order[[1]]
  } else {
    cluster_order <- NULL
  }

  seu <- readRDS(seu_path)
  seu$scna <- factor(seu$scna)
  levels(seu$scna)[levels(seu$scna) == ""] <- "diploid"

  # if(!is.null(cluster_order)){
  #
  # 	group.by = unique(cluster_order$resolution)
  #
  # 	cluster_order <-
  # 		cluster_order %>%
  # 		dplyr::mutate(order = dplyr::row_number()) %>%
  # 		dplyr::filter(!is.na(clusters)) %>%
  # 		dplyr::mutate(clusters = as.character(clusters))
  #
  # 	seu@meta.data$clusters = seu@meta.data[[group.by]]
  #
  # 	seu_meta <- seu@meta.data %>%
  # 		tibble::rownames_to_column("cell") %>%
  # 		dplyr::left_join(cluster_order, by = "clusters") %>%
  # 		dplyr::select(-clusters) %>%
  # 		dplyr::rename(clusters = phase) %>%
  # 		identity()
  #
  # 	phase_levels = phase_levels[phase_levels %in% unique(seu_meta$clusters)]
  #
  # 	seu_meta <-
  # 		seu_meta %>%
  # 		tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
  # 		dplyr::arrange(clusters, new_clusters) %>%
  # 		dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
  # 		tibble::column_to_rownames("cell") %>%
  # 		identity()
  #
  # 	seu@meta.data <- seu_meta[rownames(seu@meta.data),]
  #
  # }

  pearls_plot <- plot_distribution_of_clones_pearls(seu, seu_name = sample_id, var_x = "scna", var_y = var_y) +
    labs(title = sample_id)

  n_scnas <- n_distinct(pearls_plot$data$scna)

  plot_width <- 3 * n_scnas

  plot_path <- glue("results/{sample_id}_pearls.pdf")

  #
  ggsave(plot_path, pearls_plot, width = plot_width, height = 6)

  return(
    plot_path
  )
}

plot_clone_cc_plots <- function(seu_path, large_clone_comparisons = NULL, scna_of_interest = "1q", var_y = "phase_level", summary_function = c("total", "partial"), kept_phases = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m"), labeled_values = c("g2_m"), var_y_levels = NULL) {
	
	summary_function = match.arg(summary_function)
	
  tumor_id <- str_extract(seu_path, "SRR[0-9]*")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  seu <- readRDS(seu_path)
  
  var_y_levels <- var_y_levels %||% levels(seu@meta.data[[var_y]]) %||% levels(factor(as.numeric(seu@meta.data[[var_y]])))
  
  var_y_colors <-
  	scales::hue_pal()(length(var_y_levels)) |> 
  	set_names(var_y_levels)

  seu$scna[seu$scna == ""] <- ".diploid"
  seu$scna <- factor(seu$scna)

  if(!is.null(large_clone_comparisons)){
  	clone_comparisons <- names(large_clone_comparisons[[sample_id]])
  	clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
  	retained_clones <- clone_comparison %>%
  		str_extract("[0-9]_v_[0-9]") %>%
  		str_split("_v_", simplify = TRUE)
  	
  } else {
  	retained_clones = sort(unique(seu$clone_opt)) |> 
  		set_names()
  }
  
  seu <- seu[, seu$clone_opt %in% retained_clones]

  cc_data <- FetchData(seu, c(var_y, "G2M.Score", "S.Score", "Phase", "scna", "phase_level"))
  
  centroid_data <-
    cc_data %>%
    { if ("phase_level" %in% names(.)) dplyr::filter(., phase_level %in% kept_phases) else . } %>%
    dplyr::group_by(.data[[var_y]], scna) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score), n_cells = dplyr::n()) %>%
  	group_by(scna) %>%
    mutate(percent = proportions(n_cells) * 100) %>%
    dplyr::mutate(centroid = "centroids") %>%
  	dplyr::mutate({{var_y}} := factor(.data[[var_y]], levels = levels(cc_data[[var_y]]))) |> 
  	mutate(percent = round(percent, digits = 1)) |>
  	# dplyr::mutate(labeled = percent) |>
  	dplyr::mutate(labeled = ifelse(.data[[var_y]] %in% labeled_values, as.character(.data[[var_y]]), NA)) |>
    identity()
  
  # stacked bar plot
  bar_plot_data <-
  centroid_data |>
  	mutate({{var_y}} := fct_rev(.data[[var_y]])) |>
  	# mutate(scna = fct_rev(scna)) |>
  	identity()
  
  bar_plot <-
  bar_plot_data |> 
  	ggplot(aes(y = scna, x = percent, fill = .data[[var_y]], label = percent)) +
  	geom_bar(position = "fill", color = "black", stat = "identity", width = 0.6) +
  	# geom_text_repel(
  	# 	position = ggpp::position_fillnudge(vjust = 0.5, y = 0.4),
  	# 	direction = 'y'
  	# ) +
  	scale_y_discrete(labels = function(x) str_wrap(x, width = 20), limits = rev) +
  	# scale_fill_discrete(direction = -1) + 
  	scale_fill_manual(values = var_y_colors) + 
  	# guides(fill = "none") + 
  	NULL

  if("phase_level" %in% names(cc_data)) cc_data <- dplyr::filter(cc_data, phase_level %in% kept_phases)

  centroid_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[var_y]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[var_y]], size = percent), alpha = 1, shape = 23, colour = "black") +
  	geom_text_repel(
  		data = centroid_data,
  		color = "black",
  		size = 3,
  		aes(x = mean_x, y = mean_y, label = labeled),
  		segment.size = 0.5,
  		point.padding = 0.5,
  		force             = 0.5,
  		nudge_x           = 0.5,
  		arrow = arrow(length = unit(0.1, "npc")),
  		# arrow = arrow(),
  		min.segment.length = 0) +
  	scale_fill_manual(values = var_y_colors, drop = FALSE) + 
    # geom_text_repel(data = centroid_data, aes(x = mean_x, y = mean_y, label = .data[[var_y]]), size = 4, color = "black") +
    facet_wrap(~scna, ncol = 1) +
    guides(
    	"fill" = FALSE,
    	color = guide_legend(override.aes = list(size = 5))
    	) +
    labs(title = sample_id) +
  	# scale_size(range = c(0, 10), breaks = c(5,10,20,40,80)) +
  	scale_size_continuous(range = c(0,10), limits = c(0, 100), breaks = c(1,10,30,100)) + 
  	scale_x_continuous(
  		limits = c(-0.75, 1),
  		breaks = seq(-0.5, 1, 0.5)
  		) + 
  	scale_y_continuous(
  		limits = c(-0.5, 1.5),
  		breaks = seq(-0.5, 1.5, 0.5)
  		) + 
  	theme(
  		strip.background = element_blank(),
  		strip.text.x = element_blank(),
  		legend.title = element_text(size = 16),
  		legend.text = element_text(size = 16)
  	) +
  	NULL

  plot_list <- list("centroid" = centroid_plot, "bar" = bar_plot)
  
  mypatch <- wrap_plots(plot_list, ncol = 1) + 
  	plot_layout(heights = c(2,1)) + 
  	plot_layout(axis_titles = "collect")
  
  return(mypatch)
}

plot_cc <- function(object = NULL, group.by = "gene_snn_res.0.6", assay = "gene",
                    label = "_filtered_", color.by = "batch", mytitle = "Title",
                    faceted = TRUE, pt.size = 1, phase_levels = c(
                      "g1", "g1_s",
                      "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia",
                      "other", "s_star"
                    ), kept_phases = NULL) {
  kept_phases <- kept_phases %||% phase_levels
  heatmap_features <- object@misc$markers[[group.by]][["presto"]]

  if (is.factor(object@meta.data[[group.by]])) {
    cluster_order <- levels(object@meta.data[[group.by]]) %>%
      set_names(.)
  } else {
    cluster_order <- unique(object@meta.data[[group.by]]) %>%
      set_names(.)
  }


  object@meta.data[[group.by]] <- factor(object@meta.data[[group.by]],
    levels = cluster_order
  )
  group_by_clusters <- object@meta.data[[group.by]]
  object@meta.data$clusters <- names(cluster_order[group_by_clusters])
  object@meta.data$clusters <- factor(object@meta.data$clusters,
    levels = unique(setNames(names(cluster_order), cluster_order)[levels(object@meta.data[[group.by]])])
  )
  heatmap_features <- heatmap_features %>%
    dplyr::arrange(Cluster) %>%
    group_by(Cluster) %>%
    slice_head(n = 6) %>%
    dplyr::filter(Gene.Name %in%
      rownames(GetAssayData(object, layer = "gene", slot = "scale.data"))) %>%
    dplyr::filter(Gene.Name %in% VariableFeatures(object)) %>%
    identity()
  heatmap_features <- heatmap_features %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)
  object_heatmap <- ggplotify::as.ggplot(seu_complex_heatmap(object,
    features = heatmap_features$Gene.Name, group.by = c(
      "G2M.Score",
      "S.Score", "clusters"
    ), col_arrangement = c("clusters"),
    cluster_rows = FALSE, column_split = sort(object@meta.data$clusters),
    row_split = rev(heatmap_features$Cluster), row_title_rot = 0,
  )) + labs(title = mytitle) + theme()
  labels <- data.frame(
    clusters = unique(object[[]][["clusters"]]),
    label = unique(object[[]][["clusters"]])
  ) %>% identity()
  cc_data <- FetchData(object, c(
    "clusters", "G2M.Score", "S.Score",
    "Phase", color.by
  ))
  centroid_data <- cc_data %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()
  centroid_plot <- cc_data %>% ggplot(aes(
    x = S.Score, y = G2M.Score,
    group = .data[["clusters"]], color = .data[[color.by]]
  )) +
    geom_point(size = pt.size) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_point(
      data = centroid_data,
      aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]),
      size = 6, alpha = 0.7, shape = 23, colour = "black"
    ) +
    NULL
  facet_cell_cycle_plot <- cc_data %>% ggplot(aes(
    x = S.Score,
    y = G2M.Score, group = .data[["clusters"]], color = .data[[color.by]]
  )) +
    geom_point(size = pt.size) +
    geom_point(
      data = centroid_data,
      aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]),
      size = 6, alpha = 0.7, shape = 23, colour = "black"
    ) +
    facet_wrap(~ .data[["clusters"]], ncol = 2) +
    theme_light() +
    geom_label(
      data = labels, aes(label = label), x = max(cc_data$S.Score) +
        0.05, y = max(cc_data$G2M.Score) - 0.1, hjust = 1,
      vjust = 1, inherit.aes = FALSE
    ) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    NULL
  appender <- function(string) str_wrap(string, width = 40)
  if (faceted) {
    return(facet_cell_cycle_plot)
  } else {
    return(centroid_plot)
  }
}


montage_images <- function(plot_files, sample_id, numbat_dir = "numbat_sridhar", tile = "6", label = "montage") {
  #
  plot_images <- magick::image_read(plot_files, density = 600)

  my_montage <- magick::image_montage(plot_images, tile = tile, geometry = "800x", shadow = FALSE)

  # my_montage <- image_montage(plot_images, geometry = c('x200+10+10', 'x800+10+10', 'x100+10+10'), tile = '3x', shadow = FALSE)

  montage_path <- glue("results/{numbat_dir}/{sample_id}_{label}.pdf")

  image_write(my_montage, format = "pdf", montage_path)

  return(montage_path)
}


make_table_s09 <- make_table_s07

make_table_s10 <- make_table_s07


#' Collate per-sample summary PDFs from multiple figure targets
#'
#' @param sample_id SRR sample identifier
#' @param clone_tree_files Character vector of debranched clone tree PDF paths
#' @param fig_s03a_files Character vector of unfiltered numbat heatmap PDF paths
#' @param numbat_expression_files Character vector of numbat expression PDF paths
#' @param numbat_bulk_clone_files Character vector of numbat bulk clone PDF paths
#' @param tile magick montage tile layout (default "4x")
#' @param density Render density in DPI (default 300)
#' @return Path to the combined per-sample PDF
#' @export
collate_sample_summary <- function(sample_id,
                                   clone_tree_files,
                                   fig_s03a_files,
                                   numbat_expression_files,
                                   numbat_bulk_clone_files,
                                   tile = "4x",
                                   density = 300) {

  clone_tree_files <- unlist(clone_tree_files)
  fig_s03a_files <- unlist(fig_s03a_files)
  numbat_expression_files <- unlist(numbat_expression_files)
  numbat_bulk_clone_files <- unlist(numbat_bulk_clone_files)

  karyogram <- glue("results/{sample_id}_karyogram.pdf")
  clone_trees <- clone_tree_files[str_detect(clone_tree_files, sample_id)]
  s03a <- fig_s03a_files[str_detect(fig_s03a_files, sample_id)]
  expression <- numbat_expression_files[str_detect(numbat_expression_files, sample_id)]
  bulk_clones <- numbat_bulk_clone_files[str_detect(numbat_bulk_clone_files, sample_id)]

  # Fallback: if no fig_s03a, use unfiltered numbat heatmaps directly from disk
  if (length(s03a) == 0) {
    nb_dir <- glue("results/numbat_sridhar/{sample_id}")
    fallback <- c(
      glue("{nb_dir}/{sample_id}_unfiltered.pdf"),
      glue("{nb_dir}/{sample_id}_unfiltered_scna_var.pdf"),
      glue("{nb_dir}/{sample_id}_filtered.pdf"),
      glue("{nb_dir}/{sample_id}_filtered_scna_var.pdf")
    )
    s03a <- fallback[file.exists(fallback)]
  }

  # Build named file list for panel labels
  panel_files <- c()
  panel_labels <- c()

  if (file.exists(karyogram)) {
    panel_files <- c(panel_files, karyogram)
    panel_labels <- c(panel_labels, "Karyogram")
  }
  if (length(clone_trees) > 0) {
    panel_files <- c(panel_files, clone_trees)
    panel_labels <- c(panel_labels,
      paste0("Clone tree", if (length(clone_trees) > 1)
        paste0(" (", seq_along(clone_trees), ")") else ""))
  }
  if (length(s03a) > 0) {
    panel_files <- c(panel_files, s03a)
    s03a_labels <- basename(s03a) |>
      str_remove(paste0(sample_id, "_")) |>
      str_remove("\\.pdf$") |>
      str_replace_all("_", " ") |>
      stringr::str_to_title()
    panel_labels <- c(panel_labels, s03a_labels)
  }
  if (length(expression) > 0) {
    panel_files <- c(panel_files, expression)
    panel_labels <- c(panel_labels, "Expression")
  }
  if (length(bulk_clones) > 0) {
    panel_files <- c(panel_files, bulk_clones)
    panel_labels <- c(panel_labels, "Bulk clones")
  }

  if (length(panel_files) == 0) {
    warning("No files found for sample ", sample_id)
    return(NULL)
  }

  # Read and annotate each panel
  annotated <- list()
  for (i in seq_along(panel_files)) {
    img <- magick::image_read_pdf(panel_files[i], density = density)
    img <- magick::image_annotate(img, panel_labels[i],
                                  size = 40, color = "black",
                                  gravity = "north", weight = 700,
                                  boxcolor = "white")
    annotated[[i]] <- img
  }
  plot_images <- do.call(c, annotated)

  summary_montage <- magick::image_montage(plot_images, tile = tile,
                                           geometry = "800x800+10+10",
                                           shadow = FALSE)

  # Add sample ID as title
  summary_montage <- magick::image_annotate(summary_montage, sample_id,
                                            size = 60, color = "black",
                                            gravity = "north", weight = 700,
                                            boxcolor = "white")

  out_path <- glue("results/{sample_id}_summary.pdf")
  magick::image_write(summary_montage, format = "pdf", path = out_path)

  return(out_path)
}