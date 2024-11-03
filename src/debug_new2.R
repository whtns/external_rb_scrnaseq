
library(targets)
source("packages.R")
source("functions.R")

# fig_02 ------------------------------

debug(plot_fig_s07_08)
# # debug(plot_clone_cc_plots)
# debug(make_clone_distribution_figure)
# debug(plot_distribution_of_clones_across_clusters)

tar_load(c("integrated_seus_1q", "integrated_seus_16q", "cluster_orders", "integrated_seus_2p", "integrated_seus_6p"))

fig_s07 <- plot_fig_s07_08(integrated_seus_2p, cluster_orders, plot_path = "results/fig_s07.pdf")

fig_s08 <- plot_fig_s07_08(integrated_seus_6p, cluster_orders, plot_path = "results/fig_s08.pdf")


fig_s04 <- plot_fig_s04_06(integrated_seus_1q, cluster_orders, plot_path = "results/fig_s04.pdf")

fig_s06 <- plot_fig_s04_06(integrated_seus_16q, cluster_orders, plot_path = "results/fig_s06.pdf")


plot_paths_1q <- imap(integrated_seus_1q, ~make_clone_distribution_figure(.x, cluster_orders,
					 																					 height = 12, width = 20, plot_path = tempfile(pattern = .y, fileext = ".pdf"), heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch")))
qpdf::pdf_combine(plot_paths_1q, "results/fig_s04.pdf") |> 
	browseURL()


plot_paths_16q <- imap(integrated_seus_16q, ~make_clone_distribution_figure(.x, cluster_orders,
																																					height = 12, width = 20, plot_path = tempfile(pattern = .y, fileext = ".pdf"), heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch")))

qpdf::pdf_combine(plot_paths_16q, "results/fig_s06.pdf") |> 
	browseURL()




seu <- readRDS("output/seurat/integrated_1q/SRR13884249_filtered_seu.rds")
