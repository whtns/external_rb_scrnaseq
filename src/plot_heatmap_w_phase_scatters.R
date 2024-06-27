source("packages.R")
source("functions.R")

library(targets)
tar_load(c("filtered_seus", "cluster_orders"))

# SRR14800534 ------------------------------

test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800534_filtered_seu.rds", cluster_order = cluster_orders["SRR14800534"], mygene = "TFF1", height = 10, width = 16, equalize_scna_clones = FALSE)

browseURL("results/SRR14800534_filtered_heatmap_phase_scatter_patchwork.pdf")

# SRR13884246 ------------------------------

test0 <- plot_seu_marker_heatmap("output/seurat/SRR13884246_filtered_seu.rds", mygene = "TFF1", height = 10, width = 16, equalize_scna_clones = FALSE, group.by = "SCT_snn_res.0.4")

browseURL("results/SRR14800534_filtered_heatmap_phase_scatter_patchwork.pdf")

# filtered plot files w/ cluster names ------------------------------
filtered_plot_files <- map2(filtered_seus, cluster_orders, plot_seu_marker_heatmap, group.by = "SCT_snn_res.0.6", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered.pdf")

# filtered plot files 1.2 ------------------------------
filtered_plot_files <- map(filtered_seus, plot_seu_marker_heatmap, group.by = "SCT_snn_res.1.2", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered_1.2.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered_1.2.pdf")

# filtered plot files 1 ------------------------------
filtered_plot_files <- map(filtered_seus, plot_seu_marker_heatmap, group.by = "SCT_snn_res.1", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered_1.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered_1.pdf")

# filtered plot files 0.8 ------------------------------
filtered_plot_files <- map(filtered_seus, plot_seu_marker_heatmap, group.by = "SCT_snn_res.0.8", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered_0.8.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered_0.8.pdf")

# filtered plot files 0.6 ------------------------------
filtered_plot_files <- map(filtered_seus, plot_seu_marker_heatmap, group.by = "SCT_snn_res.0.6", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered_0.6.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered_0.6.pdf")

# filtered plot files 0.4 ------------------------------
filtered_plot_files <- map(filtered_seus, plot_seu_marker_heatmap, group.by = "SCT_snn_res.0.4", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered_0.4.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered_0.4.pdf")

# filtered plot files 0.2 ------------------------------
filtered_plot_files <- map(filtered_seus, plot_seu_marker_heatmap, group.by = "SCT_snn_res.0.2", equalize_scna_clones = FALSE)

qpdf::pdf_combine(unlist(filtered_plot_files), "results/heatmap_phase_scatter_patchworks_filtered_0.2.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_filtered_0.2.pdf")


# end integrated ------------------------------

test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800535_filtered_seu.rds", cluster_orders$SRR14800535, mygene = "TFF1", height = 10, width = 16, equalize_scna_clones = FALSE)

browseURL("results/SRR14800535_filtered_heatmap_phase_scatter_patchwork.pdf")

test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800536_wo_15q_19q_filtered_seu.rds", cluster_orders$SRR14800536, mygene = "TFF1", height = 10, width = 16)

browseURL("results/SRR14800536_filtered_heatmap_phase_scatter_patchwork.pdf")

SRR14800534_SRR14800535_SRR14800536_wo_15q_cluster_order = c("G1" = 0,
                                                             "S" = 1,
                                                             "G2" = 4,
                                                             "S*" = 5,
                                                             "post-mito" = 2,
                                                             "HSP"= 3)

# undebug(plot_seu_marker_heatmap)
test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds", cluster_order = SRR14800534_SRR14800535_SRR14800536_wo_15q_cluster_order, group.by = "integrated_snn_res.0.2", assay = "integrated", mygene = "TFF1", height = 10, width = 16)

SRR14800534_SRR14800535_SRR14800536_wo_15q_cluster_order = c("G1" = 0,
                                                             "G1" = 5,
                                                             "G1-S" = 1,
                                                             "S" = 2,
                                                             "S*" = 7,
                                                             "G2" = 4,
                                                             "post-mito" = 6,
                                                             "HSP"= 3)

# undebug(plot_seu_marker_heatmap)
test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds", cluster_order = SRR14800534_SRR14800535_SRR14800536_wo_15q_cluster_order, group.by = "integrated_snn_res.0.4", assay = "integrated", mygene = "TFF1", height = 10, width = 16)

browseURL(test0)

test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800534_SRR14800535_SRR14800536_wo_15q_seu.rds", group.by = "integrated_snn_res.0.4", assay = "integrated", mygene = "TFF1", height = 10, width = 16)

browseURL(test0)

test0 <- plot_seu_marker_heatmap("output/seurat/SRR14800536_filtered_seu.rds", cluster_orders$SRR14800536, mygene = "TFF1", height = 10, width = 16, equalized_scna_clones = FALSE)

browseURL("results/SRR14800536_filtered_heatmap_phase_scatter_patchwork.pdf")

# regressed plot files ------------------------------
regressed_plot_files <- map(regressed_seus, plot_seu_marker_heatmap, label = "_regressed_")
qpdf::pdf_combine(unlist(regressed_plot_files), "results/heatmap_phase_scatter_patchworks_regressed.pdf")

browseURL("results/heatmap_phase_scatter_patchworks_regressed.pdf")

# # combined ------------------------------
#
# # combined_meta <-
#   read_csv("results/adata_obs.csv") %>%
#   dplyr::mutate(leiden = as.factor(leiden)) %>%
#   select(c("leiden", "G2M_score", "S_score")) %>%
#   ggplot(aes(x = `S_score`, y = `G2M_score`, color = leiden)) +
#   geom_point(size = 0.5) +
#   facet_wrap(~leiden, ncol = 2) +
#   theme_light()
