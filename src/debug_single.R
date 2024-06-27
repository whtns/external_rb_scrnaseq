source("packages.R")
source("functions.R")
library(targets)

# make_numbat_heatmaps ------------------------------

# debug(make_numbat_heatmaps)
# debug(plot_phylo_heatmap)
# debug(plot_numbat)
# debug(plot_phylo_heatmap)
# 
debug(make_numbat_heatmaps)

test0 <- make_numbat_heatmaps("output/seurat/SRR27187899_filtered_seu.rds", "output/numbat_sridhar/SRR27187899_numbat.rds", p_min = 0.5, line_width = 0.1, extension = "_unfiltered")
