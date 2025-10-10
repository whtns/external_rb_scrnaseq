# Plotting and annotation functions (1)

#' Create a numbat-related plot visualization
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Parameter for myannot
#' @param mytitle Plot title
#' @param sort_by Character string (default: "scna")
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

plot_numbat <- function(nb, myseu, myannot, mytitle, sort_by = "scna", ...) {
  
  clone_levels <- sort(unique(nb$clone_post$clone_opt))
  clone_pal <- setNames(scales::hue_pal()(length(clone_levels)), clone_levels)
  scna_levels <- sort(unique(myannot$scna))
  scna_pal <- setNames(scales::hue_pal()(length(scna_levels)), scna_levels)
  myheatmap <- nb$plot_phylo_heatmap(
    pal_clone = clone_pal,
    pal_annot = scna_pal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = sort_by,
    annot_bar_width = 1,
    raster = FALSE,
    ...
  ) +
    labs(title = mytitle) +
    theme(legend.position = "none")
  myheatmap[[2]] <- myheatmap[[2]] + theme(legend.position = "none", axis.text.x = element_blank())
  return(myheatmap)
}

#' Create a numbat-related plot visualization
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Parameter for myannot
#' @param mytitle Plot title
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_numbat_w_phylo <- function(nb, myseu, myannot, mytitle, ...) {
  
  celltypes <- tibble::tibble(cell = gsub("\\.", "-", rownames(myseu@meta.data)), type = myseu@meta.data$type)
  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")
  mypal <- c("1" = "gray", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3")
  nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = TRUE,
    annot_bar_width = 1,
    ...
  ) + labs(title = mytitle)
}

safe_plot_numbat <- safely(plot_numbat, otherwise = NA_real_)
safe_plot_numbat_w_phylo <- safely(plot_numbat_w_phylo, otherwise = NA_real_)

#' Create a plot visualization
#'
#' @param phylo_plot_output Parameter for phylo plot output
#' @param chrom Character string (default: "1")
#' @param p_min Parameter for p min
#' @return ggplot2 plot object
#' @export
plot_variability_at_SCNA <- function(phylo_plot_output, chrom = "1", p_min = 0.9) {
  
  plot_input <- phylo_plot_output
  plot_input$seg <- factor(plot_input$seg, levels = str_sort(unique(plot_input$seg), numeric = TRUE))
  plot_input <- plot_input[order(plot_input$seg), ]
  plot_input$cnv_state <- dplyr::case_when(
    plot_input$cnv_state == "amp" ~ "gain",
    plot_input$cnv_state == "bamp" ~ "balanced gain",
    plot_input$cnv_state == "del" ~ "loss",
    plot_input$cnv_state == "loh" ~ "CNLoH",
    TRUE ~ plot_input$cnv_state
  )
  p_cnv_plot <- plot_input %>%
    dplyr::group_by(seg) %>%
    dplyr::filter(!all(is.na(LLR))) %>%
    ggplot(aes(x = cell_index, y = p_cnv)) +
    geom_point(aes(color = cnv_state), size = 0.1, alpha = 0.1) +
    geom_hline(aes(yintercept = p_min), linetype = 'dashed', color = "grey") +
    scale_x_reverse() +
    scale_color_manual(values = c(
      "gain" = "#7f180f",
      "balanced gain" = "pink",
      "loss" = "#010185",
      "CNLoH" = "#387229"
    )) +
    labs(color = "SCNA state", fill = "SCNA", y = "Probability of SCNA") +
    facet_wrap(~seg) +
    geom_tile(aes(y = -0.2, height = 0.1, fill = scna)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  return(p_cnv_plot)
}

#' Load or read data from file
#'
#' @param myseus Parameter for myseus
#' @return Loaded data object
#' @export
read_regress_save <- function(myseus) {
  
  seu <- myseus[[sample_id]] %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
    RunPCA(features = VariableFeatures(.), nfeatures.print = 10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters(resolution = 0.15) %>%
    RunUMAP(dims = 1:10, min.dist = 0.01)
  saveRDS(seu, myseus)
}

#' Load or read data from file
#'
#' @param myseus Parameter for myseus
#' @return Loaded data object
#' @export
read_unregress_cc_save <- function(myseus) {
  
  seu <- myseus[[sample_id]] %>%
    seuratTools::clustering_workflow(resolution = c(0.2, 0.4))
  # Uncomment below to save the object
  # saveRDS(seu, myseus)
}

read_cluster_dictionary <- function(cluster_dictionary_path = "data/cluster_dictionary.csv") {
  cluster_dictionary <- read_tsv(cluster_dictionary_path) %>%
    split(.$sample_id)

  return(cluster_dictionary)
}

assign_phase_cluster_at_resolution <- function(seu_path = NULL, cluster_order = NULL, assay = "SCT", resolution = 1, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star")) {
  #

	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
	cluster_order <- cluster_order[[file_id]]

  seu <- readRDS(seu_path)

  # start loop ------------------------------


  if (!is.null(cluster_order)) {
    single_cluster_order <- cluster_order[[resolution]]

    group.by <- unique(single_cluster_order$resolution)
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
      seu |>
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
  }
  return(saveRDS(seu, seu_path))
  # return(seu_path)
}

calculate_clone_distribution <- function(seu_path = NULL, cluster_order = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 5, width = 9, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters", pairwise = TRUE) {
  kept_phases <- kept_phases %||% phase_levels

  file_id <- fs::path_file(seu_path)
  
  tumor_id <- str_extract(seu_path, "SRR[0-9]*")
  
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  
  message(file_id)
  cluster_order <- cluster_order[[file_id]]

  full_seu <- readRDS(seu_path)

  full_seu$scna[full_seu$scna == ""] <- ".diploid"
  full_seu$scna <- factor(full_seu$scna)

  # subset by retained clones ------------------------------
  clone_comparisons <- names(large_clone_comparisons[[sample_id]])
  clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
  retained_clones <- clone_comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split("_v_", simplify = TRUE)

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  plot_path <- glue("results/numbat_sridhar/{sample_id}_clone_distribution.pdf")
  table_path <- glue("results/numbat_sridhar/{sample_id}_clone_distribution.xlsx")

  pdf(plot_path, height = height, width = width)

  pairwise_seu_tables <- list()
  for (resolution in names(cluster_order)) {
    # start loop ------------------------------

    seu <- full_seu[, full_seu$clone_opt %in% retained_clones]

    if (!is.null(cluster_order)) {
      single_cluster_order <- cluster_order[[resolution]]

      group.by <- unique(single_cluster_order$resolution)
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

    all_seu_plot <- plot_distribution_of_clones_across_clusters(seu, seu_name = sample_id, var_x = "scna", var_y = "clusters")

    if (pairwise) {
      pairwise_seu_plots <- list()

      scna_clones <- unique(sort(as.factor(seu@meta.data$scna)))

      pairwise_clone_vectors <-
        bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
        t() %>%
        as.data.frame() %>%
        as.list() %>%
        map(as.character) %>%
        identity()

      names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~ paste(., collapse = "_v_"))


      pairwise_res <- imap(pairwise_clone_vectors, make_pairwise_plots, seu, sample_id, var_y = "phase_level")

      pairwise_seu_plots <- map(pairwise_res, "plot") |>
        map(~ {
          .x + labs(title = sample_id, subtitle = resolution)
        })

      pairwise_seu_tables[resolution] <- map(pairwise_res, "table")

      print(pairwise_seu_plots)
    }
  }

  all_seu_plot <- plot_distribution_of_clones_across_clusters(seu, seu_name = sample_id, var_x = "scna", var_y = "clusters", plot_type = "clone")

  dev.off()

  table_path <- writexl::write_xlsx(pairwise_seu_tables, table_path)

  return(list("table" = table_path, "plot" = plot_path))
}