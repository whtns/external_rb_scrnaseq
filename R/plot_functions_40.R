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
plot_effect_of_filtering <- function(unfiltered_seu_path, filtered_seu_path = NULL, group.by = "gene_snn_res.0.2", cluster_dictionary,
                                     mito_threshold = c(5, 10),
                                     nCount_threshold = c(500, 1000),
                                     nFeature_threshold = c(500, 1000),
                                     low_hypoxia_seu_path = NULL,
                                     plot_path = NULL) {

  sample_id <- str_extract(unfiltered_seu_path, "SRR[0-9]*")

  if (is.null(plot_path)) {
    fs::dir_create("results/effect_of_filtering")
    plot_path <- glue::glue("results/effect_of_filtering/{sample_id}_effect_of_filtering.pdf")
  }

  unfiltered_seu <- readRDS(unfiltered_seu_path)
  original_filtered_seu <- if (!is.null(filtered_seu_path)) readRDS(filtered_seu_path) else NULL

  # Replace missing/empty scna with differential GT_opt labels
  # For each unique GT_opt value, show only the segments not in its nearest parent clone.
  # e.g. if clone A = "1q", clone B = "1q,2p", clone C = "1q,2p,13q":
  #   A -> "1q",  B -> "[1]+2p",  C -> "[2]+13q"
  differential_gt_labels <- function(gt_vals) {
    unique_gts <- unique(gt_vals[!is.na(gt_vals) & gt_vals != ""])
    if (length(unique_gts) == 0) {
      result <- gt_vals
      result[is.na(gt_vals) | gt_vals == ""] <- "diploid"
      return(result)
    }
    seg_sets <- setNames(
      lapply(unique_gts, function(gt) sort(trimws(strsplit(gt, ",")[[1]]))),
      unique_gts
    )
    sorted_gts <- names(seg_sets)[order(lengths(seg_sets))]
    # For each clone, find the closest proper-subset parent
    parent_of <- setNames(rep(NA_character_, length(sorted_gts)), sorted_gts)
    for (i in seq_along(sorted_gts)) {
      gt <- sorted_gts[i]; segs <- seg_sets[[gt]]
      best_len <- -1L
      for (j in seq_len(i - 1L)) {
        cand <- sorted_gts[j]; cand_segs <- seg_sets[[cand]]
        if (length(cand_segs) > best_len && all(cand_segs %in% segs)) {
          parent_of[gt] <- cand; best_len <- length(cand_segs)
        }
      }
    }
    clone_idx <- setNames(seq_along(sorted_gts), sorted_gts)
    display <- setNames(character(length(sorted_gts)), sorted_gts)
    for (gt in sorted_gts) {
      parent <- parent_of[gt]
      if (is.na(parent)) {
        display[gt] <- paste(seg_sets[[gt]], collapse = ",")
      } else {
        new_segs <- sort(setdiff(seg_sets[[gt]], seg_sets[[parent]]))
        display[gt] <- paste0("[", clone_idx[parent], "]+", paste(new_segs, collapse = ","))
      }
    }
    result <- gt_vals
    for (gt in sorted_gts) result[!is.na(result) & result == gt] <- display[gt]
    result[is.na(gt_vals) | (!is.na(gt_vals) & gt_vals == "")] <- "diploid"
    result
  }

  fill_scna_from_gt_opt <- function(seu) {
    if (!"GT_opt" %in% colnames(seu@meta.data)) return(seu)
    if (!"scna" %in% colnames(seu@meta.data)) {
      needs_gt <- rep(TRUE, nrow(seu@meta.data))
    } else {
      needs_gt <- is.na(seu@meta.data$scna) | seu@meta.data$scna == ""
    }
    if (!any(needs_gt)) return(seu)
    if (!"scna" %in% colnames(seu@meta.data)) seu@meta.data$scna <- NA_character_
    gt_diff <- differential_gt_labels(seu@meta.data$GT_opt[needs_gt])
    seu@meta.data$scna[needs_gt] <- gt_diff
    seu
  }
  unfiltered_seu <- fill_scna_from_gt_opt(unfiltered_seu)
  if (!is.null(original_filtered_seu)) original_filtered_seu <- fill_scna_from_gt_opt(original_filtered_seu)

  plot_list <- list()

  # use filter_reason / filter_keep columns if present, otherwise fall back to cluster_dictionary
  if (all(c("filter_reason", "filter_keep") %in% colnames(unfiltered_seu@meta.data))) {
    filtered_seu <- unfiltered_seu[, unfiltered_seu$filter_keep]
  } else {
    removed_clusters <- cluster_dictionary[[sample_id]] |>
      dplyr::filter(remove == 1) |>
      dplyr::pull(`gene_snn_res.0.2`)
    filtered_seu <- unfiltered_seu[, !unfiltered_seu$gene_snn_res.0.2 %in% removed_clusters]
  }

  # threshold sweep ------------------------------
  meta <- unfiltered_seu@meta.data
  sweep_grid <- expand.grid(
    mito_threshold = mito_threshold,
    nCount_threshold = nCount_threshold,
    nFeature_threshold = nFeature_threshold
  )
  sweep_grid$n_kept <- apply(sweep_grid, 1, function(r) {
    sum(
      meta$percent.mt < r[["mito_threshold"]] &
      meta$nCount_gene > r[["nCount_threshold"]] &
      meta$nFeature_gene > r[["nFeature_threshold"]]
    )
  })
  sweep_grid$pct_kept <- sweep_grid$n_kept / nrow(meta) * 100
  sweep_grid$mito_threshold <- factor(sweep_grid$mito_threshold)
  sweep_grid$nCount_threshold <- factor(sweep_grid$nCount_threshold)
  sweep_grid$nFeature_threshold <- factor(sweep_grid$nFeature_threshold)

  p_sweep <- ggplot(sweep_grid, aes(x = nCount_threshold, y = nFeature_threshold, fill = pct_kept)) +
    geom_tile() +
    geom_text(aes(label = n_kept), size = 2.5) +
    facet_wrap(~ mito_threshold, labeller = label_both) +
    scale_fill_viridis_c(name = "% kept", limits = c(0, 100)) +
    labs(title = glue::glue("{sample_id} — threshold sweep"), x = "nCount_threshold", y = "nFeature_threshold") +
    theme_bw()

  tmp_path <- tempfile(fileext = ".pdf")
  plot_list["threshold_sweep"] <- tmp_path
  ggsave(tmp_path, p_sweep, height = 4, width = 3 * length(mito_threshold))

  # abbreviation markers ------------------------------
  # Skip marker computation for very large objects — presto::wilcoxauc
  # (C++) crashes the process and cannot be caught by tryCatch.
  # Samples with >50k cells are too large to run safely.
  safe_plot_markers <- function(seu, grp = group.by, max_cells = 50000L, ...) {
    cat("safe_plot_markers: ncol =", ncol(seu), "\n")
    if (ncol(seu) > max_cells) {
      warning("Skipping marker plot — too many cells (", ncol(seu), ")")
      return(NULL)
    }
    # Assay5 + joined layers required by seuratTools::find_all_markers
    if (!is(seu[["gene"]], "Assay5")) {
      seu[["gene"]] <- tryCatch(as(seu[["gene"]], "Assay5"), error = function(e) seu[["gene"]])
    }
    seu[["gene"]] <- tryCatch(JoinLayers(seu[["gene"]]), error = function(e) seu[["gene"]])
    if (nrow(GetAssayData(seu, assay = "gene", layer = "data")) == 0) {
      seu <- NormalizeData(seu, assay = "gene")
    }
    tryCatch(
      seuratTools::plot_markers(seu, grp, ...),
      error = function(e) {
        warning("plot_markers failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  p_unfiltered <- safe_plot_markers(unfiltered_seu, marker_method = "wilcox", return_plotly = FALSE, hide_technical = "all", num_markers = 5)
  gc()
  p_filtered   <- safe_plot_markers(filtered_seu,   marker_method = "wilcox", return_plotly = FALSE, hide_technical = "all", num_markers = 5)
  gc()

  if (!is.null(p_unfiltered) && !is.null(p_filtered)) {
    (p_unfiltered + labs(title = "unfiltered") + theme(legend.position = "none") + NULL) +
      (p_filtered + labs(title = "filtered") + theme(axis.title.y = element_blank()) + NULL) +
      plot_annotation(title = sample_id) +
      plot_layout(guides = "collect")

    tmp_path <- tempfile(fileext = ".pdf")
    plot_list["abbreviation_markers"] <- tmp_path
    ggsave(tmp_path, height = 6, width = 9)
  }

  # filter_reason umap ------------------------------
  if ("filter_reason" %in% colnames(unfiltered_seu@meta.data)) {
    reason_levels <- c("clone_opt_na", "qc_fail", "cluster_remove", "malat1", "manual_exclude")
    unfiltered_seu@meta.data$filter_reason <- factor(unfiltered_seu@meta.data$filter_reason, levels = reason_levels)
    reason_cols <- c(scales::hue_pal()(length(reason_levels)), "grey80") |>
      set_names(c(reason_levels, "kept"))

    # sweep dimplots: re-evaluate qc_fail per threshold combo, carry forward other reasons
    sweep_plots <- lapply(seq_len(nrow(sweep_grid)), function(i) {
      r <- sweep_grid[i, ]
      qc_fail_sweep <-
        meta$percent.mt >= as.numeric(as.character(r$mito_threshold)) |
        meta$nCount_gene <= as.numeric(as.character(r$nCount_threshold)) |
        meta$nFeature_gene <= as.numeric(as.character(r$nFeature_threshold))
      sweep_reason <- dplyr::case_when(
        !is.na(unfiltered_seu@meta.data$filter_reason) & unfiltered_seu@meta.data$filter_reason != "qc_fail" ~ as.character(unfiltered_seu@meta.data$filter_reason),
        qc_fail_sweep ~ "qc_fail",
        TRUE ~ NA_character_
      )
      unfiltered_seu@meta.data$sweep_reason_label <- factor(
        tidyr::replace_na(sweep_reason, "kept"),
        levels = c(reason_levels, "kept")
      )
      DimPlot(unfiltered_seu, group.by = "sweep_reason_label") +
        scale_color_manual(values = reason_cols, limits = names(reason_cols), na.value = "grey90") +
        labs(
          title = glue::glue("mito<{r$mito_threshold}"),
          subtitle = glue::glue("nCount>{r$nCount_threshold} nFeature>{r$nFeature_threshold}"),
          colour = "filter reason"
        )
    })

    # all sweep plots on a single page: ncol = nCount thresholds, rows = mito x nFeature combos
    tmp_path <- tempfile(fileext = ".pdf")
    n_sweep_rows <- length(mito_threshold) * length(nFeature_threshold)
    sweep_ok <- tryCatch({
      pdf(tmp_path, height = 3 * n_sweep_rows, width = 4.5 * length(nCount_threshold))
      on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
      grid <- patchwork::wrap_plots(sweep_plots, ncol = length(nCount_threshold)) +
        patchwork::plot_annotation(title = sample_id) +
        patchwork::guide_area() +
        patchwork::plot_layout(guides = "collect")
      print(grid)
      dev.off()
      on.exit(NULL)
      TRUE
    }, error = function(e) {
      if (!is.null(dev.list())) dev.off()
      warning("sweep_dimplots failed: ", conditionMessage(e))
      FALSE
    })
    if (sweep_ok) plot_list["sweep_dimplots"] <- tmp_path

    reason_labels <- tidyr::replace_na(as.character(unfiltered_seu@meta.data$filter_reason), "kept")
    unfiltered_seu@meta.data$filter_reason_label <- factor(reason_labels, levels = c(reason_levels, "kept"))

    tmp_path <- tempfile(fileext = ".pdf")
    fr_ok <- tryCatch({
      pdf(tmp_path, height = 3, width = 4.5)
      on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
      print(DimPlot(unfiltered_seu, group.by = "filter_reason_label") +
        scale_color_manual(values = reason_cols, limits = names(reason_cols), na.value = "grey90") +
        labs(title = sample_id, colour = "filter reason"))
      dev.off()
      on.exit(NULL)
      TRUE
    }, error = function(e) {
      if (!is.null(dev.list())) dev.off()
      warning("filter_reason_umap failed: ", conditionMessage(e))
      FALSE
    })
    if (fr_ok) plot_list["filter_reason_umap"] <- tmp_path
  }

  # dimplot grids ------------------------------
  scna_cols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["scna"]])))

  unfiltered_seu@meta.data$scna <- vec_split_label_line(unfiltered_seu@meta.data$scna, 3)
  filtered_seu@meta.data$scna   <- vec_split_label_line(filtered_seu@meta.data$scna, 3)

  cluster_labels <- unique(as.character(unfiltered_seu@meta.data[[group.by]]))
  cluster_labels <- stringr::str_sort(cluster_labels, numeric = TRUE)
  group_cols  <- scales::hue_pal()(length(cluster_labels)) |> set_names(cluster_labels)

  n_cols <- if (!is.null(original_filtered_seu)) 3L else 2L

  # row 1: scna
  scna_plots <- list(
    DimPlot(unfiltered_seu, group.by = "scna", cols = scna_cols) + labs(title = "unfiltered"),
    DimPlot(filtered_seu,   group.by = "scna", cols = scna_cols) + labs(title = "filtered")
  )
  if (!is.null(original_filtered_seu)) {
    original_filtered_seu@meta.data$scna <- vec_split_label_line(original_filtered_seu@meta.data$scna, 3)
    scna_plots[[3]] <- DimPlot(original_filtered_seu, group.by = "scna", cols = scna_cols) + labs(title = "re-clustered")
  }

  # row 2: cluster
  cluster_plots <- list(
    DimPlot(unfiltered_seu, group.by = group.by, cols = group_cols) + labs(title = "unfiltered"),
    DimPlot(filtered_seu,   group.by = group.by, cols = group_cols) + labs(title = "filtered")
  )
  if (!is.null(original_filtered_seu)) {
    cluster_plots[[3]] <- DimPlot(original_filtered_seu, group.by = group.by, cols = group_cols) + labs(title = "re-clustered")
  }

  dimplot_grid <- patchwork::wrap_plots(c(scna_plots, cluster_plots), ncol = n_cols) +
    patchwork::plot_annotation(title = sample_id)

  tmp_path <- tempfile(fileext = ".pdf")
  plot_list["abbreviation_umaps"] <- tmp_path
  ggsave(tmp_path, dimplot_grid, height = 6, width = 4.5 * n_cols)

  # hypoxia filtering dimplots ------------------------------
  # Show which cells are dropped between filtered_seu and seus_low_hypoxia.
  # Uses original_filtered_seu (pipeline seu with UMAP) as the "before" object.
  if (!is.null(low_hypoxia_seu_path) && file.exists(low_hypoxia_seu_path)) {
    low_hypoxia_seu <- tryCatch(readRDS(low_hypoxia_seu_path), error = function(e) NULL)
    if (!is.null(low_hypoxia_seu)) {
      source_seu <- if (!is.null(original_filtered_seu)) original_filtered_seu else filtered_seu
      kept_barcodes <- colnames(low_hypoxia_seu)
      source_seu@meta.data$hypoxia_filter_label <- ifelse(
        colnames(source_seu) %in% kept_barcodes, "kept_low_hypoxia", "dropped_high_hypoxia"
      )
      hypoxia_filter_cols <- c(kept_low_hypoxia = "#2166ac", dropped_high_hypoxia = "#d73027")

      p_hypoxia_source <- DimPlot(source_seu, group.by = "hypoxia_filter_label") +
        scale_color_manual(values = hypoxia_filter_cols) +
        labs(title = glue::glue("{sample_id} filtered → hypoxia split"), colour = NULL)

      hypoxia_plot_list <- list(p_hypoxia_source)

      if ("hypoxia_score" %in% colnames(source_seu@meta.data)) {
        p_hypoxia_score <- FeaturePlot(source_seu, features = "hypoxia_score") +
          labs(title = "hypoxia score (filtered seu)")
        hypoxia_plot_list[[2]] <- p_hypoxia_score
      }

      p_low_hypoxia_scna <- DimPlot(low_hypoxia_seu, group.by = "scna") +
        labs(title = "low hypoxia seu (scna)")

      hypoxia_plot_list[[length(hypoxia_plot_list) + 1]] <- p_low_hypoxia_scna

      hypoxia_grid <- patchwork::wrap_plots(hypoxia_plot_list, ncol = length(hypoxia_plot_list)) +
        patchwork::plot_annotation(title = glue::glue("{sample_id} — hypoxia filtering"))

      tmp_path <- tempfile(fileext = ".pdf")
      hyp_ok <- tryCatch({
        ggsave(tmp_path, hypoxia_grid, height = 4, width = 4.5 * length(hypoxia_plot_list))
        TRUE
      }, error = function(e) {
        warning("hypoxia_filtering_dimplots failed: ", conditionMessage(e))
        FALSE
      })
      if (hyp_ok) plot_list["hypoxia_filtering_dimplots"] <- tmp_path
    }
  }

  # Filter to only valid, non-empty PDF files before combining
  valid_pdfs <- Filter(function(f) !is.null(f) && file.exists(f) && file.size(f) > 100L,
                       unlist(plot_list))
  if (length(valid_pdfs) == 0) {
    warning("No valid PDF pages produced for ", sample_id)
    return(NULL)
  }
  plot_path <- tryCatch(
    qpdf::pdf_combine(valid_pdfs, plot_path),
    error = function(e) {
      warning("pdf_combine failed: ", conditionMessage(e))
      NULL
    }
  )

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

  # Drop cells that were not in regressed_seu (AddMetaData gives NA for missing cells)
  cells_with_regressed <- !is.na(filtered_seu@meta.data[[colnames(regressed_meta)[1]]])
  filtered_seu <- filtered_seu[, cells_with_regressed]

  if (inherits(filtered_seu[["gene"]], "Assay5")) {
    filtered_seu <- JoinLayers(filtered_seu, assay = "gene")
  }

  valid_metavars <- colnames(regressed_meta)[
    sapply(colnames(regressed_meta), function(m) {
      dplyr::n_distinct(filtered_seu@meta.data[[m]], na.rm = TRUE) > 1
    })
  ]

  if (length(valid_metavars) == 0) {
    return(NA_character_)
  }

  filtered_seu <- tryCatch({
    find_all_markers(filtered_seu, valid_metavars)
  }, error = function(e) {
    warning("Skipping regression effect plot due to marker calculation failure: ", e$message)
    return(NULL)
  })

  if (is.null(filtered_seu)) {
    return(NA_character_)
  }

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

