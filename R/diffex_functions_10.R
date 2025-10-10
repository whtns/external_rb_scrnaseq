# Diffex Functions (10)

#' Load or read data from file
#'
#' @param mp_file File path
#' @return Loaded data object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

read_mps <- function(mp_file = "/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS") {
  mps <- readRDS(mp_file)

  mps <- purrr::map(mps, ~ purrr::set_names(.x, janitor::make_clean_names(names(.x))))
}

#' Perform differential expression analysis
#'
#' @param oncoprint_input_by_scna_unfiltered Parameter for oncoprint input by scna unfiltered
#' @return List object
#' @export
tally_num_diffex <- function(oncoprint_input_by_scna_unfiltered) {
  symbol_tally <- purrr::map_depth(oncoprint_input_by_scna_unfiltered, 2, ~ {
    table(.x$symbol)
  })
  sample_tally <- purrr::map_depth(oncoprint_input_by_scna_unfiltered, 2, ~ {
    table(.x$sample_id)
  })
  return(list("symbol" = symbol_tally, "sample" = sample_tally))
}

#' Perform clustering analysis
#'
#' @param diffex_1 Parameter for diffex 1
#' @param diffex_2 Parameter for diffex 2
#' @param new_col_name Color specification
#' @return Function result
#' @export
annotate_cluster_membership <- function(diffex_1, diffex_2, new_col_name) {
  diffex_2 <-
    diffex_2 %>%
    dplyr::ungroup() %>%
    dplyr::select(clone_comparison, symbol) %>%
    dplyr::mutate({{ new_col_name }} := TRUE) %>%
    identity()

  diffex_1 <-
    diffex_1 %>%
    dplyr::ungroup() %>%
    dplyr::left_join(diffex_2, by = c("clone_comparison", "symbol"))

  return(diffex_1)
}

#' Perform tally kooi candidates operation
#'
#' @param cis_diffex_clones Character string (default: "results/diffex_bw_clones_large_in_segment_by_chr.xlsx")
#' @param trans_diffex_clones Character string (default: "results/diffex_bw_clones_trans_by_chr.xlsx")
#' @return List object
#' @export
tally_kooi_candidates <- function(cis_diffex_clones = "results/diffex_bw_clones_large_in_segment_by_chr.xlsx", trans_diffex_clones = "results/diffex_bw_clones_trans_by_chr.xlsx") {
  #

  cis_diffex_clones <-
    cis_diffex_clones %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = cis_diffex_clones) %>%
    map(dplyr::filter, !is.na(kooi_region))

  cis_diffex_clones <-
    cis_diffex_clones[map_lgl(cis_diffex_clones, ~ (nrow(.x) > 0))] %>%
    dplyr::bind_rows(.id = "chr")

  trans_diffex_clones <-
    trans_diffex_clones %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = trans_diffex_clones) %>%
    map(dplyr::filter, !is.na(kooi_region))

  trans_diffex_clones <-
    trans_diffex_clones[map_lgl(trans_diffex_clones, ~ (nrow(.x) > 0))] %>%
    dplyr::bind_rows(.id = "chr")

  return(list("cis" = cis_diffex_clones, "trans" = trans_diffex_clones))
}

#' Perform collect study metadata operation
#'
#' @return Function result
#' @export
collect_study_metadata <- function() {
  #
  
#' Perform retrieve cell stats operation
#'
#' @param seu_path File path
#' @return Function result
#' @export
retrieve_cell_stats <- function(seu_path) {
    seu <- readRDS(seu_path)

    stats <- seu@meta.data[c("nCount_gene", "nFeature_gene", "percent.mt")] %>%
      tibble::rownames_to_column("cell")

    return(stats)
  }

  seus <-
    dir_ls("output/seurat/", regexp = "\\/SRR[0-9]*_seu.rds") %>%
    set_names(str_extract(., "SRR[0-9]*"))

  # collin ------------------------------

  collin_cell_stats <- seus[c("SRR13633759", "SRR13633760", "SRR13633761", "SRR13633762")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # field ------------------------------

  field_cell_stats <- seus[c("SRR17960480", "SRR17960481", "SRR17960482", "SRR17960483", "SRR17960484")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # wu ------------------------------

  wu_cell_stats <- seus[c("SRR13884240", "SRR13884241", "SRR13884242", "SRR13884243", "SRR13884244", "SRR13884245", "SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # yang ------------------------------

  yang_cell_stats <- seus[c("SRR14800534", "SRR14800535", "SRR14800536", "SRR14800537", "SRR14800538", "SRR14800539", "SRR14800540", "SRR14800541", "SRR14800542", "SRR14800543")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  liu_cell_stats <- seus[c("SRR27187899", "SRR27187900", "SRR27187901", "SRR27187902")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # combined ------------------------------

  study_cell_stats <- dplyr::bind_rows(list("collin" = collin_cell_stats, "field" = field_cell_stats, "wu" = wu_cell_stats, "yang" = yang_cell_stats, "liu" = liu_cell_stats), .id = "study")

  write_csv(study_cell_stats, "results/study_cell_stats.csv")

  return(study_cell_stats)
}

