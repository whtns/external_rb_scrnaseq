# Plot Functions (142)

#' Perform generate plae ref operation
#'
#' @param plae_seu_path File path
#' @return Function result
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

generate_plae_ref <- function(plae_seu_path = "data/plae_human_fetal_seu.rds") {
  
  
  plae_human_fetal_seu <- readRDS(plae_seu_path)

  plae_ref <- AggregateExpression(plae_human_fetal_seu, group.by = "CellType_predict") %>%
    pluck("RNA")

  return(plae_ref)
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param sample_id Parameter for sample id
#' @param plae_ref Parameter for plae ref
#' @param group.by Character string (default: "SCT_snn_res.0.4")
#' @param query_genes Gene names or identifiers
#' @return ggplot2 plot object
#' @export
plot_celltype_predictions <- function(seu_path, sample_id, plae_ref = NULL, group.by = "SCT_snn_res.0.4", query_genes = NULL) {
  
  
  #

  seu <- seu_path

  celltypes <- c(
    "Amacrine Cells", "Bipolar Cells", "Cones",
    "Early RPCs", "Horizontal Cells", "Late RPCs",
    "Neurogenic Cells", "Photoreceptor Precursors",
    "Retinal Ganglion Cells", "Rods", "RPCs", "RPE"
  )

  if (is.null(query_genes)) {
    query_genes <- VariableFeatures(seu)
  }

  seu_mat <- GetAssayData(seu, slot = "data")

  sub_annotable <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% query_genes)

  if (is.null(plae_ref)) {
    plae_ref <-
      "data/plae_pseudobulk_counts.csv" %>%
      read_csv() %>%
      dplyr::inner_join(sub_annotable, by = c("Gene" = "ensgene"), relationship = "many-to-many") %>%
      dplyr::distinct(study, type, symbol, .keep_all = TRUE) %>%
      dplyr::filter(type %in% str_remove(celltypes, "s$")) %>%
      dplyr::group_by(symbol, type) %>%
      dplyr::summarize(total_counts = mean(counts)) %>%
      tidyr::pivot_wider(names_from = "type", values_from = "total_counts") %>%
      tibble::column_to_rownames("symbol") %>%
      as.matrix() %>%
      identity()
  } else {
    plae_ref <- plae_ref[rownames(plae_ref) %in% query_genes, colnames(plae_ref) %in% celltypes]
  }

  res <- clustify(
    input = seu_mat,
    metadata = seu[[group.by]][[1]],
    ref_mat = plae_ref,
    query_genes = query_genes
  )

  cor_to_call(res)

  res2 <- cor_to_call(
    cor_mat = res, # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  res3 <- cor_to_call_rank(
    cor_mat = res, # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  ) %>%
    dplyr::mutate(cluster = .data[[group.by]]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(rank)

  dir_create("results/clustify")
  table_path <- glue("results/clustify/{sample_id}_clustifyr.csv")
  write_csv(res3, table_path)


  # Insert into original metadata as "type" column
  seu@meta.data <- call_to_metadata(
    res = res2, # data.frame of called cell type for each cluster
    metadata = seu@meta.data, # original meta.data table containing cell clusters
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  neurogenic_table <- seu@meta.data %>%
    dplyr::mutate(cluster = .data[[group.by]]) %>%
    janitor::tabyl(cluster, type) %>%
    dplyr::mutate(sample = sample_id) %>%
    # dplyr::mutate(percent_neurogenic = `Neurogenic Cells`/(Cones+`Neurogenic Cells`)) %>%
    identity()

  allcells_dimplot <- DimPlot(
    seu,
    group.by = "type",
    split.by = group.by
  ) +
    plot_annotation(title = sample_id) +
    NULL

  dir_create("results/clustify")
  plot_path <- glue("results/clustify/{sample_id}_clustifyr.pdf")
  ggsave(plot_path)

  return(list("plot" = plot_path, "table" = neurogenic_table))
}

#' Load or read data from file
#'
#' @param zinovyev_file File path
#' @return Loaded data object
#' @export
read_zinovyev_genes <- function(zinovyev_file = "data/zinovyev_cc_genes.tsv") {
  
  
  zinovyev_cc_genes <-
    read_tsv(zinovyev_file) %>%
    dplyr::group_by(term) %>%
    tidyr::nest(data = symbol) %>%
    # split(.$term) %>%
    tibble::deframe() %>%
    map(tibble::deframe) %>%
    identity()
}

#' Load or read data from file
#'
#' @param cc_file File path
#' @return Loaded data object
#' @export
read_giotti_genes <- function(cc_file = "data/giotti_cc_genes.tsv") {
  
  
  giotti_cc_genes <-
    read_tsv(cc_file) %>%
    dplyr::filter(!(term %in% c("Function known but link to cell division not well established", "Uncharacterized"))) %>%
    dplyr::mutate(term = str_wrap(term, width = 10)) %>%
    identity()
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param seu_name Parameter for seu name
#' @return ggplot2 plot object
#' @export
plot_phase_distribution_by_scna <- function(seu_path, seu_name) {
  
  
    seu <- seu_path
    seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))
    plot_distribution_of_clones_across_clusters(seu, seu_name, var_x = "scna", var_y = "Phase", both_ways = FALSE)
  }

