# Numbat heatmap and plotting functions (5)

#' Perform numbat-related analysis
#'
#' @param seu_path File path
#' @param numbat_rds_files File path
#' @param p_min Parameter for p min
#' @param line_width Parameter for line width
#' @param extension Character string (default: "")
#' @param midline_threshold Threshold value for filtering
#' @return ggplot2 plot object
#' @export
make_numbat_heatmaps <- function(seu_path, numbat_rds_files, p_min = 0.9, line_width = 0.1, extension = "", midline_threshold = 0.4) {
  sample_id <- str_extract(seu_path, "SRR[0-9]*")
  names(numbat_rds_files) <- str_extract(numbat_rds_files, "SRR[0-9]*")
  numbat_rds_file <- numbat_rds_files[[sample_id]]
  numbat_dir <- "numbat_sridhar"
  dir_create(glue("results/{numbat_dir}/"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))
  seu <- readRDS(seu_path)
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  retained_segs <- mynb$joint_post |> 
    dplyr::mutate(at_midline = dplyr::case_when(
      dplyr::between(p_cnv, 0.3, 0.7) ~ 1,
      .default = 0
    )) |> 
    group_by(seg) |>
    dplyr::summarise(percent_at_midline = sum(at_midline)/dplyr::n()) |>
    dplyr::filter(percent_at_midline <= midline_threshold) |>
    dplyr::arrange(desc(percent_at_midline)) |>
    dplyr::pull(seg) |>
    identity()
  mynb$joint_post <- mynb$joint_post[mynb$joint_post$seg %in% retained_segs,]
  myannot <- mynb$joint_post[, c("cell")]
  seu <- seu[, colnames(seu) %in% myannot$cell]
  myannot <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, scna) %>%
    identity()
  myannot$scna[myannot$scna == ""] <- ".diploid"
  numbat_heatmap <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)[["result"]]
  if (!is.null(numbat_heatmap) && !identical(numbat_heatmap, NA_real_)) {
    heatmap_no_phylo_path <- tempfile(fileext = ".pdf")
    ggsave(heatmap_no_phylo_path, plot = numbat_heatmap, w = 10, h = 5)
    scna_variability_plot <-
      numbat_heatmap[[3]][["data"]] |>
      dplyr::left_join(myannot, by = "cell") |>
      plot_variability_at_SCNA(p_min = p_min)
    scna_var_path <- tempfile(fileext = ".pdf")
    ggsave(scna_var_path, plot = scna_variability_plot, w = 10, h = 5)
  } else {
    heatmap_no_phylo_path <- NULL
    scna_var_path <- NULL
  }
  numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min)[["result"]]
  heatmap_w_phylo_path <- tempfile(fileext = ".pdf")
  if (!is.null(numbat_heatmap_w_phylo) && !identical(numbat_heatmap_w_phylo, NA_real_)) {
    ggsave(heatmap_w_phylo_path, plot = numbat_heatmap_w_phylo, w = 10, h = 5)
  } else {
    heatmap_w_phylo_path <- NULL
  }
  plot_path <- glue("results/{numbat_dir}/{sample_id}/{sample_id}{extension}.pdf")
  pdf_inputs <- purrr::compact(list(heatmap_no_phylo_path, scna_var_path, heatmap_w_phylo_path))
  qpdf::pdf_combine(pdf_inputs, plot_path)
  return(plot_path)
}