library(Seurat)
library(clustifyr)

celltypes = c("Amacrine Cells", "Bipolar Cells", "Cones",
              "Early RPCs", "Horizontal Cells", "Late RPCs",
              "Neurogenic Cells", "Photoreceptor Precursors",
              "Retinal Ganglion Cells", "Rods", "RPCs"
)

plae_ref <- readRDS("~/Homo_sapiens/numbat/plae_ref.rds")

plae_human_fetal_seu <- readRDS("data/plae_human_fetal_seu.rds")

generat_plae_ref <- function(plae_seu_path = "data/plae_human_fetal_seu.rds"){

  plae_human_fetal_seu <- readRDS(plae_seu_path)

  plae_ref <- AggregateExpression(plae_human_fetal_seu, group.by = "CellType_predict") %>%
    pluck("RNA")

  return(plae_ref)
}

plae_ref <- AggregateExpression(plae_human_fetal_seu, group.by = "CellType_predict") %>%
  pluck("RNA")

test0 <- readRDS("/dataVolume/storage/scEiad/human_droplet_data_fetal_seu.rds")

# seu <- readRDS("output/seurat/SRR13884242_regressed_seu.rds")

seu <- readRDS("output/seurat/SRR13884242_filtered_seu.rds")

clone_2_v_3_genes <-
  large_all_diffex_clones$large_all_diffex_clones_c2bf52b0$`3_v_2_8p+_11p+` %>%
  dplyr::arrange(avg_log2FC) %>%
  dplyr::slice_head(n = 20) %>%
  pull(symbol) %>%
  identity()


# run function ------------------------------

plot_celltype_predictions <- function(seu, sample_id, plae_ref = NULL, group.by = "SCT_snn_res.0.4", query_genes = NULL) {
  # browser()

  seu <- readRDS(seu_path)

  celltypes = c("Amacrine Cells", "Bipolar Cells", "Cones",
                "Early RPCs", "Horizontal Cells", "Late RPCs",
                "Neurogenic Cells", "Photoreceptor Precursors",
                "Retinal Ganglion Cells", "Rods", "RPCs", "RPE"
  )

  if(is.null(query_genes)){
    query_genes <- VariableFeatures(seu)
  }

  seu_mat <- GetAssayData(seu, slot = "data")

  sub_annotable <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% query_genes)

  if(is.null(plae_ref)){
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
    plae_ref = plae_ref[rownames(plae_ref) %in% query_genes,colnames(plae_ref) %in% celltypes]
  }

  res <- clustify(
    input = seu_mat,
    metadata = seu[[group.by]][[1]],
    ref_mat = plae_ref,
    query_genes = query_genes
  )

  cor_to_call(res)

  res2 <- cor_to_call(
    cor_mat = res,                  # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  res3 <- cor_to_call_rank(
    cor_mat = res,                  # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  ) %>%
    dplyr::group_by(`SCT_snn_res.0.4`) %>%
    dplyr::arrange(rank)

  # return(res3)

  # Insert into original metadata as "type" column
  seu@meta.data <- call_to_metadata(
    res = res2,                     # data.frame of called cell type for each cluster
    metadata = seu@meta.data,           # original meta.data table containing cell clusters
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  allcells_dimplot <- DimPlot(
    seu,
    group.by = "type",
    split.by = group.by
  ) +
    plot_annotation(title = sample_id) +
    NULL


  return(allcells_dimplot)

}

celltype_prediction_by_var <- function(seu_path, plae_ref, group.by = "SCT_snn_res.0.4"){

  seu <- readRDS(seu_path)

  cluster_markers <-
    table_cluster_markers(seu)[[group.by]] %>%
    group_by(Cluster) %>%
    dplyr::slice_head(n = 20) %>%
    split(.$Cluster) %>%
    map(pull, Gene.Name) %>%
    identity()

  test0 <- map(cluster_markers, ~plot_celltype_predictions(seu, plae_ref, query_genes = .x))

  return(test0)
}


# plot_celltype_predictions(seu, plae_ref, query_genes = clone_2_v_3_genes)
# rgc_genes <- c("DPYSL2", "ELAVL2", "ELAVL3",  "GNG3", "LY6H", "PCM1", "RTN1", "STMN4")

seu_path <- "output/seurat/SRR14800534_regressed_seu.rds"
seu <- readRDS(seu_path)

test4 <- plot_celltype_predictions(seu, "SRR14800534", plae_ref)

# test3 <- celltype_prediction_by_var(seu_path, plae_ref)
