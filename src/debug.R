source("packages.R")
source("functions.R")
library(targets)
tar_load("large_filter_expressions")
tar_load("cluster_dictionary")
tar_load("interesting_samples")
tar_load("large_in_segment_diffex_clones")
tar_load("large_out_of_segment_diffex_clones")
tar_load("large_clone_comparisons")

filter_diffex_for_recurrence <- function(df, num_recur = 2){
  # browser()

  test0 <-
    df %>%
    # dplyr::arrange(symbol, sample_id) %>%
    dplyr::arrange(p_val_adj) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::filter(n_distinct(sample_id) >= num_recur) %>%
    identity()

  test1 <-
    test0 %>%
    dplyr::summarize(mean_FC = mean(abs(avg_log2FC))) %>%
    dplyr::slice_max(abs(mean_FC), n =10) %>%
    dplyr::inner_join(test0, by = "symbol")

  return(test1)
}

plot_recurrence <- function(test0, mytitle){
  # browser()

  test1 <-
    test0 %>%
    dplyr::arrange(mean_FC) %>%
    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) %>%
    dplyr::select(sample_id, symbol, neg_log_p_val_adj, abs_log2FC, mean_FC) %>%
    identity()

  ggplot(test1, aes(x=sample_id, y=symbol)) +
    geom_point(aes(color=neg_log_p_val_adj, size = abs_log2FC)) +
    # geom_tile(fill = NA, color = "black", linewidth = 0.5) +
    labs(title = mytitle, color = "-log \n p_adj") +
    theme_bw() +
    scale_size_continuous(
      limits = c(0.1, 0.9)
    ) +
    scale_color_continuous(
      limits = c(1, 150),
      oob = squish
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
      # legend.position="none"
    )
}

# plot drivers ------------------------------

samples_1q <- c("SRX11133594", "SRX11133592", "SRX11133593", "SRX10264526",
                "SRX14116944")

samples_2p <- c("SRX10264526", "SRX14116944", "SRX14116947", "SRX10264524")

samples_6p <- c("SRX14116944", "SRX10264524", "SRX14116947")

samples_16q <- c("SRX11133588", "SRX10264519", "SRX11133593", "SRX11133594", "SRX11133585", "SRX11133592", "SRX10264520")

# samples_16q <- c("SRX11133588", "SRX11133585", "SRX11133592")

plot_markers_featureplot <- function(sample_ids, features){
  # browser()

  seus <- glue("output/seurat/{sample_ids}_filtered_seu.rds") %>%
    set_names(str_extract(., "SRR[0-9]*")) %>%
    map(readRDS) %>%
    identity()

  markerplots <- map(seus, FeaturePlot, features = features, combine = TRUE)

  markerplots0 <- imap(markerplots, ~{.x + labs(subtitle = .y)})

  vlnplots <- map(seus, VlnPlot, features = features, group.by = "GT_opt", combine = TRUE, pt.size = 0)

  vlnplots0 <- imap(vlnplots, ~{.x +
      labs(title = .y) +
      theme(legend.position="none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank())})

  # vln_wrapped = wrap_plots(vlnplots0, ncol = 2) + plot_annotation(title = features)

  return(vlnplots0)
}

vlns_1q <- plot_markers_featureplot(samples_1q, "CENPF")

vlns_2p <- plot_markers_featureplot(samples_2p, "SOX11")

vlns_6p <- plot_markers_featureplot(samples_6p, "DEK")

vlns_16q <- plot_markers_featureplot(samples_16q, "CDT1")


vlns_1q +
  wrap_plots(nrow = 1) + plot_annotation(title = "CENPF")
ggsave("results/vln_plots_from_oncoprint_1q.pdf", width = 6, height = 8)

vlns_2p +
  wrap_plots(nrow = 1) + plot_annotation(title = "SOX11")
ggsave("results/vln_plots_from_oncoprint_2p.pdf", width = 6, height = 6)
vlns_6p
ggsave("results/vln_plots_from_oncoprint_6p.pdf", width = 6, height = 6)
vlns_16q
ggsave("results/vln_plots_from_oncoprint_16q.pdf", width = 6, height = 6)
