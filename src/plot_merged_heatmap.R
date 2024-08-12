#!/usr/bin/env Rscript
library(targets)
source("packages.R")
source("functions.R")

tar_load(c("debranched_seus_16q", "debranched_seus_1q"))

# debug("plot_merged_heatmap")
# debug(seu_gene_heatmap)
# debug("seu_complex_heatmap2")

# test0 <- plot_merged_heatmap(seu_path = "output/seurat/merged_16q_filtered_seu.rds", child_seu_paths = debranched_seus_16q) 

test0 <- plot_merged_heatmap(seu_path = "output/seurat/merged_1q_filtered_seu.rds", child_seu_paths = debranched_seus_1q) 
