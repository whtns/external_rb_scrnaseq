# Plot Functions (124)

#' Create a plot visualization
#'
#' @param mymarkers Parameter for mymarkers
#' @param seu_paths File path
#' @param plot_type Parameter for plot type
#' @param group_by Character string (default: "gene_snn_res.0.2")
#' @param extension Character string (default: "_filtered")
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

plot_putative_marker_across_samples <- function(mymarkers, seu_paths, plot_type = FeaturePlot, group_by = "gene_snn_res.0.2", extension = "_filtered") {
  print(mymarkers)

  sample_ids <- str_extract(seu_paths, "SRR[0-9]*")

  myplots <- map(seu_paths, plot_markers_in_sample, mymarkers = mymarkers, plot_type = plot_type, group_by = group_by) %>%
    set_names(sample_ids)


  myplots0 <-
    myplots %>%
    transpose() %>%
    imap(~ {
      patchwork::wrap_plots(.x) +
        plot_annotation(
          title = .y
        )
    })

  if (identical(plot_type, VlnPlot)) {
    plot_type_label <- "VlnPlot"
  } else if (identical(plot_type, FeaturePlot)) {
    plot_type_label <- "FeaturePlot"
  }


  plot_paths <- glue("results/numbat_sridhar/gene_plots/{names(myplots0)}_{plot_type_label}_{group_by}_{extension}.pdf")

  map2(plot_paths, myplots0, ~ ggsave(.x, .y, height = 8, width = 14))

  plot_path <- qpdf::pdf_combine(plot_paths, glue("results/numbat_sridhar/gene_plots/{plot_type_label}_{group_by}_{extension}.pdf"))

  fs::file_delete(plot_paths)

  return(plot_path)
}

#' Perform differential expression analysis
#'
#' @param numbat_rds_file File path
#' @param cluster_dictionary Cluster information
#' @param ident.1 Cell identities or groups
#' @param ident.2 Cell identities or groups
#' @return ggplot2 plot object
#' @export
find_diffex_bw_clusters_for_each_clone <- function(numbat_rds_file, cluster_dictionary, ident.1 = "G2M", ident.2 = "cone") {
  #
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")


  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[, !is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  
#' Perform clustering analysis
#'
#' @param clone_for_diffex Parameter for clone for diffex
#' @param seu Seurat object
#' @return Function result
#' @export
cluster_diff_per_clone <- function(clone_for_diffex, seu) {
    #
    seu0 <- seu[, seu$clone_opt == clone_for_diffex]

    Idents(seu0) <- seu0$abbreviation

    diffex <- FindAllMarkers(seu0, group.by = "abbreviation")
  }

  myclones <- sort(unique(seu$clone_opt)) %>%
    set_names(.)

  possible_cluster_diff_per_clone <- possibly(cluster_diff_per_clone)

  diffex <- map(myclones, possible_cluster_diff_per_clone, seu)

  diffex0 <- map(diffex, compact) %>%
    map(tibble::rownames_to_column, "symbol") %>%
    dplyr::bind_rows(.id = "clone") %>%
    dplyr::select(-symbol) %>%
    dplyr::rename(symbol = gene) %>%
    dplyr::mutate(sample_id = sample_id) %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::select(symbol, description, everything()) %>%
    dplyr::distinct(symbol, clone, .keep_all = TRUE) %>%
    dplyr::arrange(cluster, clone, p_val_adj) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    identity()

  diffex_path <- glue("results/{numbat_dir}/{sample_id}_clone_diffex.csv")
  write_csv(diffex0, diffex_path)


  
#' Perform differential expression analysis
#'
#' @param diffex_bw_clusters_for_each_clone Cluster information
#' @return Differential expression results
#' @export
compare_diffex_cluster_by_clone <- function(diffex_bw_clusters_for_each_clone) {
    test0 <-
      diffex_bw_clusters_for_each_clone %>%
      dplyr::group_by(clone, cluster) %>%
      dplyr::slice_head() %>%
      dplyr::arrange(cluster, clone)
  }

  diffex1 <- compare_diffex_cluster_by_clone(diffex0)

  cluster_clone <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, clone_opt, abbreviation) %>%
    tidyr::unite(cluster_clone, abbreviation, clone_opt) %>%
    dplyr::select(cell, cluster_clone) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- AddMetaData(seu, cluster_clone)

  diffex_cluster_by_clone_dotplot <-
    DotPlot(seu, features = unique(diffex1$symbol), group.by = "cluster_clone") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    coord_flip() +
    labs(title = sample_id)

  dotplot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_clone_diffex_dotplot.pdf")
  ggsave(dotplot_path, diffex_cluster_by_clone_dotplot, width = 10, height = 12)

  trend_genes <- diffex1 %>%
    slice_head(n = 1)

  diffex_cluster_by_clone_trendplot <-
    plot_gene_clone_trend(seu, trend_genes$symbol) +
    labs(title = sample_id)

  trendplot_path <- glue("results/{numbat_dir}/{sample_id}_gene_trend.pdf")
  ggsave(trendplot_path, diffex_cluster_by_clone_trendplot, width = 10, height = 30)

  return(list("diffex" = diffex_path, "plot" = dotplot_path))
}
#' Perform clustering analysis
#'
#' @param cluster_for_diffex Cluster information
#' @param seu Seurat object
#' @param group.by Parameter for group.by
#' @return Function result
#' @export
clone_diff_per_cluster <- function(cluster_for_diffex, seu, group.by) {
    #
    seu0 <- seu[, seu[["clusters"]] == cluster_for_diffex]

    Idents(seu0) <- seu0$clone_opt

    # diffex <- FindAllMarkers(seu0) %>%
    #   dplyr::rename(clone_opt = cluster)

    diffex <- imap(clone_comparisons, make_clone_comparison_integrated, seu0, mynbs, location = location)

    return(diffex)
  }

plot_figure_collage <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters") {
	kept_phases <- kept_phases %||% phase_levels
	
	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
	cluster_order <- cluster_order[[file_id]]
	
	full_seu <- readRDS(seu_path)
	
	# subset by retained clones ------------------------------
	clone_comparisons <- names(large_clone_comparisons[[sample_id]])
	clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
	retained_clones <- clone_comparison %>%
		str_extract("[0-9]_v_[0-9]") %>%
		str_split("_v_", simplify = TRUE)
	
	nb_paths <- nb_paths %>%
		set_names(str_extract(., "SRR[0-9]*"))
	
	nb_path <- nb_paths[[tumor_id]]
	
	plot_paths <- vector(mode = "list", length = length(cluster_order))
	names(plot_paths) <- names(cluster_order)
	
	file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	plot_path <- glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork.pdf")
	
	pdf(plot_path, height = height, width = width)
	
	for (resolution in names(cluster_order)) {
		# start loop ------------------------------
		
		seu <- full_seu[, full_seu$clone_opt %in% retained_clones]
		
		if (!is.null(cluster_order)) {
			single_cluster_order <- cluster_order[[resolution]]
			
			group.by <- unique(single_cluster_order$resolution)
		}
		
		if (equalize_scna_clones) {
			seu_meta <- seu@meta.data %>%
				tibble::rownames_to_column("cell")
			
			clones <- table(seu_meta$scna)
			
			min_clone_num <- clones[which.min(clones)]
			
			selected_cells <-
				seu_meta %>%
				dplyr::group_by(scna) %>%
				slice_sample(n = min_clone_num) %>%
				pull(cell)
			
			seu <- seu[, selected_cells]
		}
		
		
		if (!is.null(single_cluster_order)) {
			single_cluster_order <-
				single_cluster_order |>
				dplyr::mutate(order = dplyr::row_number()) %>%
				dplyr::filter(!is.na(clusters)) %>%
				dplyr::mutate(clusters = as.character(clusters))
			
			group.by <- unique(single_cluster_order$resolution)
			
			seu@meta.data$clusters <- seu@meta.data[[group.by]]
			
			seu_meta <- seu@meta.data %>%
				tibble::rownames_to_column("cell") %>%
				dplyr::select(-any_of(c("phase_level", "order"))) %>%
				dplyr::left_join(single_cluster_order, by = "clusters") %>%
				dplyr::select(-clusters) %>%
				dplyr::rename(phase_level = phase) %>%
				identity()
			
			phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$phase_level)]
			
			seu_meta <-
				seu_meta %>%
				tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
				dplyr::arrange(phase_level, order) %>%
				dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
				tibble::column_to_rownames("cell") %>%
				identity()
			
			seu@meta.data <- seu_meta[rownames(seu@meta.data), ]
			
			seu <-
				seu[, seu$phase_level %in% kept_phases] %>%
				find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
				identity()
			
			seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)
			
			# mysec ------------------------------
			
			heatmap_features <-
				seu@misc$markers[["clusters"]][["presto"]] %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu))
			
			tidy_eval_arrange <- function(.data, ...) {
				.data %>%
					arrange(...)
			}
			#
			single_cluster_order_vec <-
				seu@meta.data %>%
				dplyr::select(clusters, !!group.by) %>%
				dplyr::arrange(clusters, !!sym(group.by)) %>%
				dplyr::select(clusters, !!group.by) |>
				dplyr::distinct(.data[[group.by]], .keep_all = TRUE) |>
				dplyr::mutate(!!group.by := as.character(.data[[group.by]])) |>
				tibble::deframe() |>
				identity()
			
			heatmap_features[["Cluster"]] <-
				factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))
			
			heatmap_features <-
				heatmap_features %>%
				dplyr::group_by(Gene.Name) |>
				dplyr::slice_max(order_by = Average.Log.Fold.Change, n = 1) |>
				dplyr::ungroup() |>
				dplyr::arrange(Cluster, desc(Average.Log.Fold.Change)) |>
				group_by(Cluster) %>%
				slice_head(n = 5) %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
				dplyr::ungroup() %>%
				dplyr::distinct(Gene.Name, .keep_all = TRUE) |>
				identity()
		} else {
			heatmap_features <-
				seu@misc$markers[[group.by]][["presto"]]
			
			single_cluster_order <- levels(seu@meta.data[[group.by]]) %>%
				set_names(.)
			
			seu@meta.data[[group.by]] <-
				factor(seu@meta.data[[group.by]], levels = single_cluster_order)
			
			group_by_clusters <- seu@meta.data[[group.by]]
			
			seu@meta.data$clusters <- names(single_cluster_order[group_by_clusters])
			
			seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(single_cluster_order), single_cluster_order)[levels(seu@meta.data[[group.by]])]))
			
			heatmap_features <-
				heatmap_features %>%
				dplyr::arrange(Cluster) %>%
				group_by(Cluster) %>%
				slice_head(n = 6) %>%
				dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
				identity()
		}
		
		large_enough_clusters <-
			seu@meta.data %>%
			dplyr::group_by(clusters) %>%
			dplyr::count() |>
			dplyr::filter(n >= min_cells_per_cluster) %>%
			dplyr::pull(clusters)
		
		seu <- seu[, seu$clusters %in% large_enough_clusters]
		
		seu$scna[seu$scna == ""] <- ".diploid"
		seu$scna <- factor(seu$scna)
		# levels(seu$scna)[1] <- "none"
		
		giotti_genes <- read_giotti_genes()

		heatmap_features <-
			heatmap_features %>%
			dplyr::ungroup() %>%
			left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
			# select(Gene.Name, term) %>%
			dplyr::mutate(term = replace_na(term, "")) %>%
			dplyr::distinct(Gene.Name, .keep_all = TRUE)
		
		row_ha <- ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))
		
		if (!is.null(split_columns)) {
			column_split <- sort(seu@meta.data[[split_columns]])
			column_title <- unique(column_split)
		} else {
			column_split <- split_columns
			column_title <- NULL
		}
		
		seu_heatmap <- ggplotify::as.ggplot(
			seu_complex_heatmap(seu,
													features = heatmap_features$Gene.Name,
													group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
													col_arrangement = c("clusters", "scna"),
													cluster_rows = FALSE,
													column_split = column_split,
													row_split = rev(heatmap_features$Cluster),
													row_title_rot = 0,
													column_title = column_title,
													column_title_rot = 90
			)
		) +
			labs(title = sample_id) +
			theme()
		
		labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
			# dplyr::rename({{group.by}} := cluster) %>%
			identity()
		
		cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))
		
		centroid_data <-
			cc_data %>%
			dplyr::group_by(clusters) %>%
			dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
			dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
			dplyr::mutate(centroid = "centroids") %>%
			identity()
		
		centroid_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
			geom_point(size = 0.1) +
			theme_light() +
			theme(
				strip.background = element_blank(),
				strip.text.x = element_blank()
			) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
			guides(fill = "none", color = "none") +
			NULL
		
		
		facet_cell_cycle_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
			geom_point(size = 0.1) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
			facet_wrap(~ .data[["clusters"]], ncol = 2) +
			theme_light() +
			geom_label(
				data = labels,
				aes(label = label),
				# x = Inf,
				# y = -Inf,
				x = max(cc_data$S.Score) + 0.05,
				y = max(cc_data$G2M.Score) - 0.1,
				hjust = 1,
				vjust = 1,
				inherit.aes = FALSE
			) +
			theme(
				strip.background = element_blank(),
				strip.text.x = element_blank()
			) +
			# guides(color = "none") +
			NULL
		
		appender <- function(string) str_wrap(string, width = 40)
		
		labels <- data.frame(scna = unique(seu$scna), label = str_replace(unique(seu$scna), "^$", "diploid"))
		
		#
		
		clone_ratio <- janitor::tabyl(as.character(seu$scna))$percent[[2]]
		
		comparison_scna <-
			janitor::tabyl(as.character(seu$scna))[2, 1]
		
		clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
			seu,
			seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
		) + 
			coord_flip()
		
		# umap_plots <- DimPlot(full_seu, group.by = c("scna", "clusters"), combine = FALSE) %>%
		# 	# map(~(.x + theme(legend.position = "bottom"))) %>%
		# 	wrap_plots(ncol = 1)
		# full_seu$clusters
		# full_seu[[group.by]] <-
		full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
		levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
		umap_plots <- make_faded_umap_plots(full_seu, retained_clones, group_by = group.by)
		
		if (!is.null(nb_path)) {
			clone_tree_plot <-
				plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = FALSE)
			
			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = facet_cell_cycle_plot,
				"C" = clone_tree_plot,
				"D" = centroid_plot,
				"E" = plot_spacer(),
				"F" = plot_spacer(),
				"G" = plot_spacer(),
				"H" = clone_distribution_plot
			)
			
			# 
			
			layout <- "
              AAAAAACDEEGG
              AAAAAABBEEGG
              AAAAAABBFFGG
              AAAAAABBFFGG
              AAAAAAHHFFGG
              "
			
			plot_collage <- wrap_plots(collage_plots) +
				# plot_layout(widths = c(16, 4)) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		} else {
			layout <- "
              AAAAAAAAAABBBBCCCC
              AAAAAAAAAABBBBCCCC
              AAAAAAAAAABBBBDDDD
              AAAAAAAAAABBBBDDDD
              AAAAAAAAAABBBBDDDD
      "
			
			collage_plots <- list(
				"seu_heatmap" = seu_heatmap,
				"facet_cell_cycle_plot" = facet_cell_cycle_plot,
				"centroid_plot" = centroid_plot,
				"clone_distribution_plot" = clone_distribution_plot
			)
			
			plot_collage <- wrap_plots(collage_plots) +
				# plot_layout(widths = c(16, 4)) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		}
		
		print(plot_collage)
		# end loop------------------------------
	}
	
	dev.off()
	
	return(plot_path)
}

plot_fig_09_10 <- function(seu_subset, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments, plot_path = "result/dotplot_recurre.pdf", widths = rep(8, 3), heights = rep(12, 3), common_seus = NULL, ...){
	
	common_seus <- common_seus %||% seu_subset
	
	names(corresponding_seus) <- fs::path_file(unlist(corresponding_seus))
	names(corresponding_clusters_diffex) <- names(corresponding_seus)
	names(corresponding_clusters_enrichments) <- names(corresponding_seus)
	
	seu_subset <- corresponding_seus[corresponding_seus %in% seu_subset]
	
	corresponding_clusters_diffex <- corresponding_clusters_diffex[names(seu_subset)]
	
	corresponding_clusters_enrichments <- corresponding_clusters_enrichments[names(seu_subset)]
	
	enrichment_files <- compare_corresponding_enrichments(corresponding_clusters_enrichments, common_seus, plot_path = tempfile(fileext = ".pdf"), width = widths[[1]]*2, height = heights[[1]])
	
	all_dotplot <- corresponding_clusters_diffex |> 
		map_depth(2, "all") |> 
		dotplot_recurrent_genes(...)
	
	cis_dotplot <- corresponding_clusters_diffex |> 
		map_depth(2, "cis") |> 
		dotplot_recurrent_genes(...)
	
	trans_dotplot <- corresponding_clusters_diffex |> 
		map_depth(2, "trans") |> 
		dotplot_recurrent_genes(...)
	
	fig_09_10_panels <- list(
		"all" = all_dotplot,
		"cis" = cis_dotplot,
		"trans" = trans_dotplot
	) |> 
		imap(~{.x + labs(title = .y)}) %>%
		identity()
	
	tmpplots <- 
		list(
		"plot" = fig_09_10_panels, 
		"height" = heights, 
		"width" = widths) |> 
		pmap(~ggsave(tempfile(fileext = ".pdf"), plot = ..1, height = ..2, width = ..3))
	
	qpdf::pdf_combine(c(tmpplots, enrichment_files), plot_path)
	
	return(plot_path)
	
}