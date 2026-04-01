# Numbat utility functions (10)
#' Perform numbat-related analysis
#'
#' @param numbat_dir Parameter for numbat dir
#' @param kept_samples Parameter for kept samples
#' @return Function result
#' @export
retrieve_numbat_rds_files <- function(numbat_dir, kept_samples = NULL) {
  numbat_rds_files <- fs::dir_ls(numbat_dir, regexp = ".*SRR[0-9]*_numbat.rds", recurse = TRUE) %>%
    set_names(str_extract(., "SRR[0-9]*"))
  if (!is.null(kept_samples)) {
    numbat_rds_files <- numbat_rds_files[names(numbat_rds_files) %in% kept_samples]
  }
  return(sort(numbat_rds_files))
}

#' Perform retrieve seus operation
#'
#' @param seu_dir Parameter for seu dir
#' @param kept_samples Parameter for kept samples
#' @return Modified Seurat object
#' @export
retrieve_seus <- function(seu_dir, kept_samples = NULL) {
  seus <- fs::dir_ls(seu_dir, regexp = "./SRR[0-9]*_seu.rds", recurse = TRUE) %>%
    set_names(str_extract(., "SRR[0-9]*"))
  if (!is.null(kept_samples)) {
    seus <- seus[names(seus) %in% kept_samples]
  }
  return(sort(seus))
}

#' Create a numbat-related plot visualization
#'
#' @param numbat_plots Parameter for numbat plots
#' @param plot_type Character string (default: "exp_roll_clust.png")
#' @return ggplot2 plot object
#' @export
retrieve_numbat_plot_type <- function(numbat_plots, plot_type = "exp_roll_clust.png") {
  retrieved_plot_types <- map(numbat_plots, ~ set_names(.x, fs::path_file(.x))) %>%
    map(str_detect, plot_type)
  retrieved_plots <- map2(numbat_plots, retrieved_plot_types, ~ {
    .x[.y]
  }) %>%
    unlist()
}

#' Perform reroute done to results pdf operation
#'
#' @param numbat_rds_file File path
#' @param label Character string (default: "")
#' @return Function result
#' @export
reroute_done_to_results_pdf <- function(numbat_rds_file, label = "") {
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")
  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]
  numbat_pdfs <- dir_ls(glue("results/{numbat_dir}/{sample_id}/"), glob = "*.pdf")
  results_file <- glue("results/{numbat_dir}/{sample_id}{label}.pdf")
  qpdf::pdf_combine(numbat_pdfs, results_file)
  return(results_file)
}