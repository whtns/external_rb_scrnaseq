# Parameter and metadata functions (11)

#' Perform retrieve snakemake params operation
#'
#' @param numbat_rds_file File path
#' @return List object
#' @export
retrieve_snakemake_params <- function(numbat_rds_file) {
  str_extract(numbat_rds_file, "SRR[0-9]*")
  log_file <- fs::path(path_dir(numbat_rds_file), "log.txt")
  log <- read_lines(log_file)[3:26] %>%
    str_split(" = ") %>%
    transpose() %>%
    identity()
  params <- log[[2]] %>%
    set_names(log[[1]])
  return(list(sample_id, params))
}

#' Perform retrieve current param operation
#'
#' @param current_params Parameter for current params
#' @param myparam Parameter for myparam
#' @return Function result
#' @export
retrieve_current_param <- function(current_params, myparam) {
  sample_ids <- map(current_params, 1)
  param_values <- map(current_params, 2) %>%
    map(myparam) %>%
    set_names(sample_ids)
}