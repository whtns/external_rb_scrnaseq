#!/usr/bin/env Rscript

source("packages.R")
source("functions.R")

library(targets)

# debug(simplify_gt)

tar_load(cluster_dictionary)
tar_load(large_clone_simplifications)
tar_load(large_filter_expressions)

debug(plot_numbat)
debug(numbat:::plot_phylo_heatmap)

# make_numbat_plot_files("output/numbat_sridhar/SRR13884242/done.txt",
#                        cluster_dictionary,
#                        filter_expressions = large_filter_expressions,
#                        clone_simplifications = large_clone_simplifications)

make_numbat_heatmaps("output/numbat_sridhar/SRR13884242/done.txt", p_min = 0.5, line_width = 0.1, large_filter_expressions, cluster_dictionary, extension = "_filtered")
