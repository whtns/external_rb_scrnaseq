library(targets)
suppressPackageStartupMessages(source("packages.R"))
source("functions.R")
library(plotgardener)

source("src/heatmap_functions.R")

tar_load(c("cluster_orders"))

seu6 <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6.rds") |>
    add_hypoxia_score() |>
    identity()

seu <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall.rds") |>
    add_hypoxia_score()

# saveRDS(seu, "output/seurat/integrated_1q_16q/integrated_seu_1q_afterall.rds")

seu_hypoxia_0_5 <- seu6[,seu6$hypoxia_score < 0.5]

saveRDS(seu_hypoxia_0_5, "output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_hypoxia_0_5.rds")

seu_hypoxia_0_5 <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_hypoxia_0_5.rds")

seu_hypoxia_0_4 <- seu6[,seu6$hypoxia_score < 0.4]

saveRDS(seu_hypoxia_0_4, "output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_hypoxia_0_4.rds")

seu_hypoxia_0_4 <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_hypoxia_0_4.rds")

thresh_plot_path = "results/hypoxia_threshold6.pdf"
pdf(thresh_plot_path, w = 6, h = 6)
plot_hypoxia_score(seu6) + 
    geom_hline(yintercept = 0.4, color = "red") + 
    geom_hline(yintercept = 0.5, color = "red")

plot_hypoxia_score(seu_hypoxia_0_5) + 
    geom_hline(yintercept = 0.5, color = "red") + 
    ylim(0,1)

plot_hypoxia_score(seu_hypoxia_0_4) + 
    geom_hline(yintercept = 0.4, color = "red") + 
    ylim(0,1)
dev.off()
browseURL(thresh_plot_path)

# make_clone_distribution_figure ------------------------------

test1 <- make_clone_distribution_figure("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6.rds", cluster_order = cluster_orders, height = 12, width = 16, file_id = "integrated_seu_1q_afterall6.rds", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score", "batch"))

browseURL(test1)

test2 <- make_clone_distribution_figure("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_less_0_4.rds", cluster_order = cluster_orders, height = 12, width = 16, file_id = "integrated_seu_1q_afterall6_less_0_4.rds", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score"))

browseURL(test2)

test3 <- make_clone_distribution_figure("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_less_0_5.rds", cluster_order = cluster_orders, height = 12, width = 16, file_id = "integrated_seu_1q_afterall6_less_0_5.rds", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score"))

browseURL(test3)

test4 <- make_clone_distribution_figure("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_more_0_4.rds", cluster_order = cluster_orders, height = 12, width = 16, file_id = "integrated_seu_1q_afterall6_more_0_4.rds", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score"))

browseURL(test4)

test4 <- make_clone_distribution_figure("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_more_0_5.rds", cluster_order = cluster_orders, height = 12, width = 16, file_id = "integrated_seu_1q_afterall6_more_0_5.rds", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score"))

browseURL(test4)

# heatmap_marker_genes_debug ------------------------------

debug(heatmap_marker_genes_debug)

test0 <- heatmap_marker_genes_debug(seu, list("hypoxia" = mt_genes, "MT" = hypoxia_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "hypoxia2", "hypoxia1"), col_arrangement = c("hypoxia_score")) |> 
    browseURL()


heatmap_marker_genes_debug(seu[,seu$hypoxia_score < 0.5], list("hypoxia" = mt_genes, "MT" = hypoxia_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "hypoxia2", "hypoxia1"), col_arrangement = c("hypoxia_score")) |> 
    browseURL()

ggplotify::as.ggplot(
    seu_complex_heatmap(seu,
                        features = heatmap_features$Gene.Name,
                        group.by = heatmap_groups,
                        col_arrangement = c("clusters", "MT1"),
                        # col_arrangement = c("cluster", "scna"),
                        cluster_rows = FALSE,
                        column_split = column_split,
                        row_split = rev(heatmap_features$Cluster),
                        row_title_rot = 0,
                        column_title = column_title,
                        column_title_rot = 90
    )
) +
    labs(title = sample_id) +
    theme()

ggsave("~/tmp.pdf") |> 
    browseURL()

test0 <- seu@meta.data[c("hypoxia_score", "clusters")] |>
    # dplyr::arrange(clusters, MT1) |>
    dplyr::arrange(clusters, hypoxia_score) |>
    tibble::rownames_to_column("cell") |> 
    dplyr::mutate(cell = factor(cell, levels = cell)) |> 
    ggplot(aes(x = cell, y = hypoxia_score)) +
    geom_point()

test0

ggsave("~/tmp.pdf") |> 
    browseURL()
