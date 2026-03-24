# Plotting and annotation functions (4)

#' Filter data based on specified criteria
#'
#' @param numbat_rds_file File path
#' @param seus Parameter for seus
#' @param cluster_dictionary Cluster information
#' @param large_clone_simplifications Parameter for large clone simplifications
#' @param filter_expressions Parameter for filter expressions
#' @param cells_to_remove Cell identifiers or information
#' @param extension Character string (default: "")
#' @param leiden_cluster_file File path
#' @return Modified Seurat object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible
# - multiple_joins: Combine multiple joins into single join operation where possible

filter_cluster_save_seu <- function(numbat_rds_file, seus, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove, extension = "", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv") {
  
  output_plots <- list()
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")
  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]
  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))
  names(seus) <- str_extract(seus, "SRR[0-9]*")
  seu <- readRDS(seus[[sample_id]])
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  nb_clone_post <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")
  seu <- Seurat::AddMetaData(seu, nb_clone_post)
  seu <- seu[, !is.na(seu$clone_opt)]
  all_cells_meta <- seu@meta.data
  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")
  seu <- AddMetaData(seu, test0)
  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")
  low_numbat_prob_cells <- map(filter_expressions[[sample_id]], pull_cells_matching_expression, phylo_heatmap_data) %>%
    unlist()
  large_clone_simplifications <-
    tibble::enframe(large_clone_simplifications[[sample_id]], "scna", "seg") %>%
    tidyr::unnest(seg)

  scna_metadata <-
    nb_clone_post %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::rowwise() %>%
    # dplyr::mutate(scna = simplify_gt_col(GT_opt, large_clone_simplifications)) %>%
    dplyr::mutate(scna = list(simplify_gt_col(GT_opt, large_clone_simplifications))) %>%
    dplyr::mutate(scna = list(wrap_scna_labels(scna))) %>%
    tidyr::unnest(scna)  |> 
    dplyr::distinct(cell, .keep_all = TRUE)  |> 
    tibble::column_to_rownames("cell") %>%
    dplyr::mutate(scna = as.character(scna))  |> 
    identity()

    scna_metadata <- scna_metadata[colnames(seu),]

    seu <- Seurat::AddMetaData(seu, scna_metadata)


  scna_meta <- seu@meta.data
  seu <- filter_sample_qc(seu,
    mito_threshold = 10, nCount_threshold = 1000, nFeature_threshold = 1000
  )
  qc_meta <- seu@meta.data
  clusters_to_remove <-
    cluster_dictionary[[sample_id]] %>%
    dplyr::filter(remove == "1") %>%
    dplyr::pull(`gene_snn_res.0.2`)
  seu <- seu[, !(seu$gene_snn_res.0.2 %in% clusters_to_remove)]
  mysample <- sample_id
  if ("MALAT1" %in% unique(seu$abbreviation)) {
    seu <- seu[, seu$abbreviation != "MALAT1"]
  }
  seu <- seu[, !colnames(seu) %in% cells_to_remove[[sample_id]][["cell"]]]
  seu <- SCTransform(seu, assay = "gene", verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- seurat_cluster(
    seu = seu, resolution = seq(0.2, 1.0, by = 0.2),
    reduction = "pca"
  )
  seu <- find_all_markers(seu, seurat_assay = "SCT")
  filtered_seu_path <- glue("output/seurat/{sample_id}_filtered_seu2.rds")
  Project(seu) <- sample_id
  add_hash_metadata(seu = seu, filepath = filtered_seu_path)

  cell_type_meta <- seu@meta.data
  plot_filtering_timeline(all_cells_meta, scna_meta, qc_meta, cell_type_meta, sample_id)
  ggsave(glue("results/{sample_id}_filtering_timeline_new.pdf"), width = 8, height = 4)
  return(filtered_seu_path)
}

#' Filter data based on specified criteria
#'
#' @param numbat_rds_file File path
#' @param cluster_dictionary Cluster information
#' @param large_clone_simplifications Parameter for large clone simplifications
#' @param filter_expressions Parameter for filter expressions
#' @param extension Character string (default: "")
#' @return Modified Seurat object
#' @export
prep_unfiltered_seu <- function(numbat_rds_file, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, extension = "") {
  
  
  output_plots <- list()
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")
  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]
  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))
  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  nb_clone_post <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")
  seu <- Seurat::AddMetaData(seu, nb_clone_post)
  seu <- seu[, !is.na(seu$clone_opt)]
  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")
  seu <- AddMetaData(seu, test0)
  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")
  large_clone_simplifications <-
    tibble::enframe(large_clone_simplifications[[sample_id]], "scna", "seg") %>%
    tidyr::unnest(seg)
  scna_metadata <-
    nb_clone_post %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(scna = simplify_gt_col(GT_opt, large_clone_simplifications)) %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    tibble::column_to_rownames("cell") %>%
    identity()
  seu <- Seurat::AddMetaData(seu, scna_metadata)
  mysample <- sample_id
  seu <- SCTransform(seu, assay = "gene", verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- seurat_cluster(
    seu = seu, resolution = c(0.2, 0.4, 0.6),
    reduction = "pca"
  )
  seu <- find_all_markers(seu, seurat_assay = "SCT")
  unfiltered_seu_path <- glue("output/seurat/{sample_id}_unfiltered_seu.rds")

  add_hash_metadata(seu = seu, filepath = unfiltered_seu_path)
  return(unfiltered_seu_path)
}

#' Filter data based on specified criteria
#'
#' @param filtered_seu_path File path
#' @return Modified Seurat object
#' @export
regress_filtered_seu <- function(filtered_seu_path) {
  
  
  regressed_seu_path <- str_replace(filtered_seu_path, "_filtered", "_regressed")
  filtered_seu <- readRDS(filtered_seu_path)
  regressed_seu <- filtered_seu
  regressed_seu <- PercentageFeatureSet(regressed_seu, pattern = "^MT-", col.name = "percent.mt")
  regressed_seu <- SCTransform(regressed_seu, assay = "gene", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)
  regressed_seu <- RunPCA(regressed_seu, verbose = FALSE)
  regressed_seu <- RunUMAP(regressed_seu, dims = 1:30, verbose = FALSE)
  regressed_seu <- FindNeighbors(regressed_seu, dims = 1:30, verbose = FALSE)
  regressed_seu <- seurat_cluster(
    seu = regressed_seu, resolution = seq(0.2, 1.0, by = 0.2),
    reduction = "pca"
  )
  regressed_seu <- find_all_markers(regressed_seu, seurat_assay = "SCT")
  add_hash_metadata(seu = regressed_seu, filepath = regressed_seu_path)
  return(regressed_seu_path)
}


filter_sample_qc <- function(seu, mito_threshold = 5, nCount_threshold = 500, nFeature_threshold = 500) {
  #
  seu <-
    seu %>%
    subset(subset = percent.mt < mito_threshold) %>%
    subset(subset = nFeature_gene > nFeature_threshold) %>%
    subset(subset = nCount_gene > nCount_threshold) %>%
    identity()
}

#' Read all hypoxia Seurat RDS files and pull hypoxia_score metadata
#'
#' @param dir Directory to search for hypoxia Seurat files. Defaults to "output/seurat".
#' @param pattern Regex pattern to match hypoxia Seurat files. Defaults to "_seu_hypoxia.rds$".
#' @param recursive Whether to search recursively. Defaults to TRUE.
#' @return A tibble with columns: sample_id, file, cell, hypoxia_score
#' @export
read_all_hypoxia_scores <- function(files) {

    res <- purrr::map_dfr(files, function(f) {
  tryCatch({
    sample_id <- stringr::str_replace(basename(f), "^([^_]+)_.*$", "\\1")
    seu <- readRDS(f)                       # errors here are caught
    md <- seu@meta.data                     # errors here are also caught
    if (!"hypoxia_score" %in% colnames(md)) {
      warning("No hypoxia_score in meta.data for ", f)
      return(tibble::tibble(sample_id = sample_id, file = f, cell = rownames(md), hypoxia_score = NA_real_))
    }
    tibble::tibble(
      sample_id = sample_id,
      file = f,
      cell = rownames(md),
      hypoxia_score = md[["hypoxia_score"]]
    )
  }, error = function(e) {
    warning("Failed processing ", f, ": ", e$message)
    # return an empty tibble row (or whatever shape you prefer)
    tibble::tibble(sample_id = stringr::str_replace(basename(f), "^([^_]+)_.*$", "\\1"),
                   file = f,
                   cell = NA_character_,
                   hypoxia_score = NA_real_)
  })
})

    return(res)
}
