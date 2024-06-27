source("packages.R")
source("functions.R")

library(speckle)
library(limma)
library(ggraph)

library(targets)
tar_load(c("filtered_seus", "regressed_seus"))

# seu <- readRDS("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds")
#
# seu <- seu[,seu$scna != "16q-,1q+,15q+,19q-"]
#
# saveRDS(seu, "output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds")

filtered_seus <-
  c(filtered_seus, "output/seurat/SRR14800536_wo_15q_19q_filtered_seu.rds") %>%
  set_names(str_extract(., "SRR[0-9]*"))

run_speckle <- function(seu_path = "output/seurat/SRR14800534_filtered_seu.rds", group.by = "SCT_snn_res.0.6") {
  # browser()
  seu <- readRDS(seu_path)
  sample_id = str_extract(seu_path, "SRR[0-9]*")

  if(!"sample_id" %in% colnames(seu@meta.data)){
    seu <- AddMetaData(seu, sample(2, ncol(seu), replace = TRUE), "sample_id")
  }

  seu_meta <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(clusters = .data[[group.by]], group = scna, sample_id, `orig.ident`) %>%
    # dplyr::mutate(group = as.numeric(factor(group))) %>%
    dplyr::mutate(sample_group = paste0(sample_id, "_", group)) %>%
    identity()

  # Run propeller testing for cell type proportion differences between the two
  # groups
  mytable <- propeller(clusters = seu_meta$clusters, sample = seu_meta$sample_group,
            group = seu_meta$group) %>%
    tibble::rownames_to_column("cluster")

  seu_meta$clusters <- factor(seu_meta$clusters, levels = mytable$cluster)

  cell_prop_plot <- plotCellTypeProps(sample=seu_meta$clusters, clusters=seu_meta$group) +
    labs(title = sample_id)

  return(list("table" = mytable, "plot" = cell_prop_plot))
}

safe_run_speckle <- purrr::safely(run_speckle)

# run_speckle("output/seurat/SRR14800534_filtered_seu.rds")

seu_paths <- c("output/seurat/SRR14800534_filtered_seu.rds",
               "output/seurat/SRR14800535_filtered_seu.rds",
               "output/seurat/SRR14800536_wo_15q_19q_filtered_seu.rds")

# propeller results ------------------------------


run_speckle_for_set <- function(seu_paths, group.by = "SCT_snn_res.0.6"){
  # browser()
  dir_create("results/speckle")
  test0 <- map(seu_paths, safe_run_speckle, group.by = group.by)

  plotlist <- compact(map(test0, c("result", "plot")))

  pdf(glue("results/speckle/{group.by}.pdf"))
  print(plotlist)
  dev.off()

  tables <- compact(map(test0, c("result", "table")))

  write_xlsx(tables, glue("results/speckle/{group.by}.xlsx"))

}

make_clustree_for_speckle_set <- function(seu_path, assay = "SCT", resolutions = seq(0.2, 2.0, by = 0.2)){

  # browser()
  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)

  scna_clones <- unique(seu$scna)

  resolutions <-
    glue("{assay}_snn_res.{resolutions}") %>%
    set_names(.)

  speckle_output <- map(resolutions, ~safe_run_speckle(seu_path, .x)) %>%
    map("result")

  collate_clone_prop <- function(mydf){
    mydf0 <-
      mydf %>%
      dplyr::select(cluster, starts_with("PropMean")) %>%
      identity()
  }

  # speckle_proportions <- map(speckle_output, "table") %>%
    # map(dplyr::mutate, significance = ifelse(FDR < 5e-2, "significant", "not")) %>%
    # map(dplyr::select, clusters = cluster, starts_with("PropMean"), significance) %>%
    # imap(~dplyr::mutate(.x, clusters = paste0(.y, "C", clusters))) %>%
    # dplyr::bind_rows() %>%
    # identity()

  speckle_proportions <-
    map(speckle_output, "plot") %>%
    map("data") %>%
    map(dplyr::mutate, Clusters = na_if(Clusters, "")) %>%
    map(dplyr::mutate, Clusters = replace_na(Clusters, "diploid")) %>%
    map(~tidyr::pivot_wider(.x, names_from = "Clusters", values_from = "Proportions")) %>%
    map(janitor::clean_names) %>%
    # map(dplyr::select, clusters = cluster, starts_with("PropMean"), significance) %>%
    imap(~dplyr::mutate(.x, clusters = paste0(.y, "C", samples))) %>%
    dplyr::bind_rows() %>%
    identity()


  clustree_plot <- clustree::clustree(seu, assay = assay)

  clustree_plot$data <- dplyr::left_join(clustree_plot$data, speckle_proportions, by = c("node" = "clusters"))

  # clustree_plot$data <- dplyr::left_join(clustree_plot$data, speckle_signif, by = c("node" = "clusters"))

  clustree_plot$layers[[2]] <- NULL

  min_scale = round(min(speckle_proportions$x16q_1q), digits =2)

  my_breaks <- c(min_scale, 0.05, 0.1, 0.2)

  clustree_plot +
    ggraph::geom_node_point(aes(colour = x16q_1q, size = size)) +
    ggraph::geom_node_text(aes(label = cluster)) +
    labs(title = sample_id) +
    # scale_color_gradient(trans = "log", breaks = my_breaks, labels = my_breaks) +
    scale_color_gradient() +
    # scale_colour_manual(c("blue", "red")) +
    NULL

}

safe_make_clustree_for_speckle_set <- safely(make_clustree_for_speckle_set)

# test0 <- map(filtered_seus, safe_make_clustree_for_speckle_set)

undebug(make_clustree_for_speckle_set)

test2 <- safe_make_clustree_for_speckle_set("output/seurat/SRR14800534_filtered_seu.rds", assay = "SCT", resolutions = c(0.2, 0.4))

test1 <- safe_make_clustree_for_speckle_set("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds", assay = "integrated")

test1$result + labs(title = "integrated_1q+_16q-")

ggsave("results/integrated_1q+_16q-_clustree.pdf", width = 8, height = 10)

browseURL("results/integrated_1q+_16q-_clustree.pdf")




seu <- readRDS("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds")

resolutions = glue("integrated_snn_res.{seq(0.2, 2.0, by = 0.2)}")

pdf("results/cluster_and_phase.pdf")
map(resolutions, ~{
  DimPlot(seu, group.by = .x) + DimPlot(seu, group.by = "Phase")
})
dev.off()

browseURL("results/cluster_and_phase.pdf")


map(paste0("SCT_snn_res.", seq(0.2, 2.0, by = 0.2)), ~run_speckle_for_set(filtered_seus, .x))

plot_cluster_and_phase <- function(seu_path, cluster_col){
  seu <- readRDS(seu_path)
  sample_id <- str_extract(seu_path, "SRR[0-9]*")
  DimPlot(seu, group.by = c(cluster_col, "Phase")) +
    plot_annotation(title = sample_id)
}

map(filtered_seus[1:2], ~map(cluster_vars, function(x){plot_cluster_and_phase(.x, x)}))


test1 <- run_speckle("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds",
                     group.by = "integrated_snn_res.0.2")

VlnPlot(seu, group.by = "integrated_snn_res.0.4", features = c("TOP2A", "RRM2"))

excel_proportions <- map(test0, "table") %>%
  map(tibble::rownames_to_column, "cluster") %>%
  write_xlsx("results/propeller_proportions.xlsx")

browseURL("results/propeller_proportions.xlsx")
