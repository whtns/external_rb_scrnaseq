source("packages.R")
source("functions.R")

library(speckle)
library(limma)
library(ggraph)

library(targets)
tar_load(c("filtered_seus"))

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

make_clustree_for_speckle_set <- function(seu_path, mylabel = "sample_id", assay = "SCT", resolutions = seq(0.2, 2.0, by = 0.2)){

  # browser()
  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)

  scna_clones <- unique(seu$scna)

  resolutions <-
    glue("{assay}_snn_res.{resolutions}") %>%
    set_names(.)

  seu_meta <- seu@meta.data %>%
    dplyr::mutate(scna = na_if(scna, "")) %>%
    dplyr::mutate(scna = replace_na(scna, "diploid")) %>%
    identity()

  speckle_proportions <- map(resolutions, ~seu_meta[,c(.x, "scna")]) %>%
    map(set_names, c("cluster", "scna")) %>%
    map(janitor::tabyl, cluster, scna) %>%
    map(janitor::adorn_percentages) %>%
    map(dplyr::rename, samples = cluster) %>%
    imap(~dplyr::mutate(.x, clusters = paste0(.y, "C", samples))) %>%
    dplyr::bind_rows() %>%
    janitor::clean_names() %>%
    identity()


  clustree_plot <- clustree::clustree(seu, assay = assay)

  clustree_plot$data <- dplyr::left_join(clustree_plot$data, speckle_proportions, by = c("node" = "clusters"))

  clustree_plot$layers[[2]] <- NULL

  color_clustree_by_clone <- function(clustree_plot, mycolor, myclone, mylabel = "asdf"){

    myplot <- clustree_plot +
      ggraph::geom_node_point(aes(colour = .data[[myclone]], size = size)) +
      ggraph::geom_node_text(aes(label = cluster)) +
      labs(title = mylabel, subtitle = myclone) +
      scale_color_distiller(palette = mycolor, direction = 1) +
      NULL

    return(myplot)
  }

  clone_names <- colnames(speckle_proportions)[!colnames(speckle_proportions) %in% c("samples", "clusters")]

  brewer_palettes <- c("Reds", "Greens", "Blues", "Purples", "Oranges")

  clone_colors <- brewer_palettes[1:length(clone_names)] %>%
    set_names(clone_names)

  clustree_plots <- imap(clone_colors, ~color_clustree_by_clone(clustree_plot, .x, .y, mylabel = mylabel))

  pdf(glue("results/{mylabel}.pdf"), width = 8, height = 10)
  clustree_plots
  dev.off()

  return(clustree_plots)

}

safe_make_clustree_for_speckle_set <- safely(make_clustree_for_speckle_set)

debug(make_clustree_for_speckle_set)
combined_clustree <- safe_make_clustree_for_speckle_set("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds", assay = "integrated", mylabel = "integrated_1q+_16q-")

# combined_clustree$result + labs(title = "integrated_1q+_16q-")

pdf("results/integrated_1q+_16q-_clustree.pdf", width = 8, height = 10)
combined_clustree
dev.off()

combined_clustrees <- imap(filtered_seus, ~make_clustree_for_speckle_set(.x, assay = "SCT", mylabel = .y))

safe_make_clustree_for_speckle_set(filtered_seus[[1]], assay = "SCT", mylabel = names(filtered_seus)[[1]])


pdf("results/clustrees.pdf", width = 8, height = 10)
combined_clustrees
dev.off()

browseURL("results/integrated_1q+_16q-_clustree.pdf")

# heatmap ------------------------------

cluster_orders <- read_tsv("data/rb_tumor_cluster_identities_per_cell_cycle.tsv") %>%
  split(.$sample_id) %>%
  map(dplyr::select, cluster_id, `SCT_snn_res.0.6`,) %>%
  map(tibble::deframe) %>%
  identity()

speckle_seu_marker_heatmap <- function(seu, cluster_order = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", mygene =  "EZH2", label = "_filtered_", height = 14, width = 22, equalize_scna_clones = FALSE, heatmap_title = NULL) {
  # browser()

  heatmap_title <- glue("{heatmap_title} {group.by}")

  if(equalize_scna_clones){
    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell")

    clones <- table(seu_meta$scna)

    min_clone_num <- clones[which.min(clones)]

    selected_cells <-
      seu_meta %>%
      dplyr::group_by(scna) %>%
      slice_sample(n = min_clone_num) %>%
      pull(cell)

    seu <- seu[,selected_cells]

  }

  heatmap_features  <-
    table_cluster_markers(seu, assay = assay)

  if(!is.null(cluster_order)){
    # browser()

    seu@meta.data[[group.by]] <-
      factor(seu@meta.data[[group.by]], levels = cluster_order)

    group_by_clusters <- seu@meta.data[[group.by]]

    seu@meta.data$clusters <- names(cluster_order[group_by_clusters])

    seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(cluster_order), cluster_order)[levels(seu@meta.data[[group.by]])]))

    heatmap_features[[group.by]][["Cluster"]] <-
      factor(heatmap_features[[group.by]][["Cluster"]], levels = cluster_order)

    heatmap_features[[group.by]] <-
      heatmap_features[[group.by]] %>%
      dplyr::arrange(Cluster) %>%
      identity()
  } else {
    seu@meta.data$clusters <- seu@meta.data[[group.by]]
  }

  heatmap_features <-
    heatmap_features %>%
    pluck(group.by) %>%
    group_by(Cluster) %>%
    slice_head(n=6) %>%
    dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
    dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
    identity()

  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  heatmap_features <-
    heatmap_features %>%
    dplyr::ungroup() %>%
    left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
    select(Gene.Name, term) %>%
    dplyr::mutate(term = replace_na(term, "")) %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)

  row_ha = ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

  seu_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(seu,
                        features = heatmap_features$Gene.Name,
                        group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
                        col_arrangement = c("clusters", "scna"),
                        cluster_rows = FALSE
    )) +
    labs(title = heatmap_title) +
    theme()
  # return(seu_heatmap)

  # browser()
  labels <- data.frame(cluster=unique(seu[[]][["clusters"]]), label =unique(seu[[]][["clusters"]])) %>%
    # dplyr::rename({{group.by}} := cluster) %>%
    dplyr::rename(clusters = cluster) %>%
    identity()


  group_by_plot <-
    FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", mygene)) %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, color = .data[["clusters"]])) +
    geom_point(size = 0.5) +
    facet_wrap(~.data[["clusters"]], ncol = 2) +
    theme_light() +
    geom_label(data = labels, aes(label = label),
               x = Inf, y = -Inf, hjust=1, vjust=0,
               inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    # guides(color = "none") +
    NULL

  # ggplot(diamonds, aes(x, y)) + xlim(4, 10) + ylim(4, 10) +
  #   stat_bin2d(aes(fill = after_stat(density))) +
  #   scale_fill_gradient(name = "Percent",
  #                       labels = scales::percent)

  appender <- function(string) str_wrap(string, width = 40)

  labels <- data.frame(scna=unique(seu$scna), label=str_replace(unique(seu$scna), "^$", "diploid"))

  # asdf
  scna_plot <-
    FetchData(seu, c("scna", "G2M.Score", "S.Score", "Phase", mygene)) %>%
    # ggplot(aes(x = `S.Score`, y = `G2M.Score`, color = .data[["Phase"]])) +
    ggplot(aes(x = `S.Score`, y = `G2M.Score`)) +
    stat_bin2d(bins = 120, aes(fill = after_stat(density))) +
    # geom_density_2d_filled() +
    # geom_density_2d_filled() +
    facet_wrap(~.data[["scna"]], ncol = 1, labeller = as_labeller(appender)) +
    # guides(color = "none") +
    theme_light() +
    geom_label(data = labels, aes(label = label),
               x = Inf, y = -Inf, hjust=1, vjust=0,
               inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    guides(color = "none") +
    scale_fill_gradient(name = "Percent", labels = scales::percent, high = "#132B43", low = "#56B1F7", trans = "log") +
    # scale_fill_gradient(name = "Percent", labels = scales::percent, high = "#132B43", low = "#56B1F7") +
    NULL

  # scatter_plots <- FeatureScatter(seu, "S.Score", "G2M.Score", group.by = group.by) /
  #   FeatureScatter(seu, "S.Score", "G2M.Score", group.by = "Phase") /
  #   FeatureScatter(seu, "S.Score", "G2M.Score", group.by = "scna")

  layout <- "
            AAAAAAAAAAAABBBBB###
            AAAAAAAAAAAABBBBBCCC
            AAAAAAAAAAAABBBBBCCC
            AAAAAAAAAAAABBBBBCCC
            AAAAAAAAAAAABBBBB###
            "

  wrap_plots(seu_heatmap, group_by_plot, scna_plot) +
    # plot_layout(widths = c(16, 4)) +
    plot_layout(design = layout) +
    NULL

  ggsave(glue("results/{label}heatmap_phase_scatter_patchwork.pdf"), height = height, width = width)

  # seu_meta <- seu@meta.data %>%
  #   # dplyr::mutate(leiden = factor(leiden)) %>%
  #   identity()
  #
  # plot_groups <- levels(seu[[]][[group.by]]) %>%
  #   set_names(.)
  #
  # pal_cats <- seu_meta[["scna"]] %>%
  #   unique()
  #
  # mypal <-
  #   scales::hue_pal()(length(pal_cats)) %>%
  #   rev() %>%
  #     set_names(pal_cats)
  #
  # distplots <-
  #   map(plot_groups, ~{
  #   seu_meta %>%
  #     dplyr::filter(`SCT_snn_res.0.6` == .x) %>%
  #     ggplot(aes(x = .data[["SCT_snn_res.0.6"]], fill = .data[["scna"]])) +
  #     geom_bar(position = "fill") +
  #     labs(title = .x) +
  #     scale_fill_manual(values = mypal) +
  #     theme(
  #       axis.title.x = element_blank(),
  #       axis.title.y = element_blank(),
  #       axis.text.x = element_blank(),
  #       legend.position = "none"
  #     )
  # })

  # plot_distribution_of_clones_across_clusters(seu, sample_id, "SCT_snn_res.0.6", "scna", both_ways = FALSE)
  #
  # ggsave(glue("results/{sample_id}{label}scna_proportions.pdf"), height = 4, width = 6)

}


# SRR14800534_SRR14800535_SRR14800536_wo_15q_cluster_order = c("G1" = 0,
#                                                              "G1" = 1,
#                                                              "S" = 2,
#                                                              "G2" = 3,
#                                                              "HSP"= 4,
#                                                              "S*" = 5)

# undebug(speckle_seu_marker_heatmap)
# test0 <- speckle_seu_marker_heatmap("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds", cluster_order = SRR14800534_SRR14800535_SRR14800536_wo_15q_cluster_order, group.by = "integrated_snn_res.0.2", assay = "integrated", mygene = "TFF1", height = 10, width = 16)


seu <- readRDS(glue("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds"))
heatmap_plots <- list()
for (i in glue("integrated_snn_res.{seq(0.2, 2.0, by = 0.2)}")){
  print(i)

  heatmap_plots[i] <- speckle_seu_marker_heatmap(seu, group.by = i, assay = "integrated", mygene = "TFF1", height = 10, width = 16, heatmap_title = "integrated_1q+_16q-", label = i)

}

qpdf::pdf_combine(heatmap_plots, "results/heatmap_patchworks.pdf")

browseURL("results/heatmap_patchworks.pdf")

