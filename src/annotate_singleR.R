library(targets)
suppressPackageStartupMessages(source("packages.R"))
source("functions.R")
library(plotgardener)


seu_path <- "output/seurat/integrated_1q/integrated_seu_1q_complete.rds"
seu <- readRDS(seu_path)

library(scran)

debug(annotate_cell_cycle_singleR)

seu <- annotate_cell_cycle_singleR(seu)

# seu$seurat_Phase <- seu$Phase
# seu$Phase <- seu$singleR_phase

saveRDS(seu, seu_path)

tar_load(cluster_orders)

fig_2a_c <- make_clone_distribution_figure(seu_path, cluster_orders,
																					 height = 12, width = 20, plot_path = "results/fig_29a_c.pdf", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))


plot_distribution_of_clones_across_clusters(
	seu2,
	seu_name = glue("test0"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
)
