plot_hypoxia_score <- function(seu) {
    seu@meta.data |> 
        tibble::rownames_to_column("cell") |> 
        dplyr::arrange(hypoxia_score) |> 
        dplyr::mutate(cell = factor(cell, levels = cell)) |> 
        ggplot(aes(x = cell, y = hypoxia_score)) +
        geom_point()
}

add_hypoxia_score <- function(seu) {
    # mt_genes <- str_subset(rownames(seu), "MT-.*")
    mt_genes <- c("MT-CO3", "MT-ND3", "MT-CYB", "MT-ATP6", "MT-CO2")
    
    # hypoxia_genes <-
    #     dplyr::filter(seu@misc$markers$clusters$presto, Cluster == "hypoxia_2") |> 
    #     slice_head(n = 100) |>
    #     dplyr::pull(Gene.Name) |>
    #     identity()
    
    hypoxia_genes <- c("BNIP3", "GAS5")
    
    nbin <- floor(length(VariableFeatures(seu)) / 100)
    
    seu <- Seurat::AddModuleScore(seu, features = list("hypoxia" = hypoxia_genes, "MT" = mt_genes), name = "hypoxia", nbin = nbin, ctrl = 100)
    
    seu$hypoxia <- seu$hypoxia1
    seu$MT <- seu$hypoxia2
    
    seu$MT <- seu$MT*-1
    
    seu$hypoxia_score <-
        rowMeans(seu@meta.data[c("hypoxia", "MT")])
    
    seu$hypoxia_score = scales::rescale(seu$hypoxia_score, c(0,1))
    
    return(seu)
}

subset_to_1q <- function(seu, file_id = NULL, slug="", ...) {
  clone_selection <- tribble(
      ~batch, ~clone_opt,
      "SRR13884249", c(1,2),
      "SRR14800534", c(2,3),
      "SRR14800535", c(2,3),
      "SRR14800536", c(2,3),
      "SRR13884246", c(1,2),
      "SRR17960484", c(1,2)
  ) |>
      tidyr::unnest(clone_opt)
  
  selected_cells <- seu@meta.data |>
      tibble::rownames_to_column("cell") |>
      dplyr::inner_join(clone_selection, by = c("batch", "clone_opt"))
  
  seu_1q <- seu[,selected_cells$cell]
  
  seu_1q$scna <-
      factor(ifelse(str_detect(seu_1q$GT_opt, pattern = "1[a-z]"), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))
  
  seu_1q <- assign_phase_clusters(seu_1q, file_id, ...)
  
  return(seu_1q)
}

select_1q_clones <- function(seu, slug="", ...) {
    
    seu_1q <- subset_to_1q(seu, ...)
    
    rds_path <-  glue("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6{slug}.rds")
    
    saveRDS(seu_1q, rds_path)
    
    # seu_1q  |> 
    #     SplitObject(split.by = "batch") |> 
    #     # imap(~saveRDS(.x, glue("output/seurat/integrated_1q_16q/{.y}_integrated_1q_after6{slug}_filtered_seu.rds"))) |> 
    #     identity()
    
    pdf_path <- make_clone_distribution_figure_debug(seu_1q, rds_path, group.bys = "clusters", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score", "batch"), ...)
    
    return(pdf_path)
    
}

subset_to_16q <- function(seu, file_id = NULL, slug="", ...) {
    clone_selection <- tribble(
        ~batch, ~clone_opt,
        "SRR14800534", c(1,2),
        "SRR14800535", c(1,2),
        "SRR14800536", c(1,2)
    ) |>
        tidyr::unnest(clone_opt)
    
    selected_cells <- seu@meta.data |>
        tibble::rownames_to_column("cell") |>
        dplyr::inner_join(clone_selection, by = c("batch", "clone_opt"))
    
    seu_16q <- seu[,selected_cells$cell]
    
    seu_16q$scna <-
        factor(ifelse(str_detect(seu_16q$GT_opt, pattern = "16[a-z]"), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))
    
    seu_16q <- assign_phase_clusters(seu_16q, file_id, ...)
    
    return(seu_16q)
}

select_16q_clones <- function(seu, file_id = NULL, slug="", ...) {
    
    seu_16q <- subset_to_16q(seu, ...)
    
    rds_path <-  glue("output/seurat/integrated_1q_16q/integrated_seu_16q_afterall6{slug}.rds")
    
    saveRDS(seu_16q, rds_path)
    
    # seu_1q  |> 
    #     SplitObject(split.by = "batch") |> 
    #     # imap(~saveRDS(.x, glue("output/seurat/integrated_1q_16q/{.y}_integrated_1q_after6{slug}_filtered_seu.rds"))) |> 
    #     identity()
    
    pdf_path <- make_clone_distribution_figure_debug(seu_16q, rds_path, group.bys = "clusters", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score", "batch"), file_id = file_id, ...)
    
    return(pdf_path)
    
}

check_integrated_cluster_numbers <- function(seu_path){
    
    seu <- readRDS(seu_path)
    
    n_clusters <- seu@meta.data[c("integrated_snn_res.0.4", "integrated_snn_res.0.6", "integrated_snn_res.0.8")] |> 
        map(dplyr::n_distinct)
    
    return(n_clusters)
}

