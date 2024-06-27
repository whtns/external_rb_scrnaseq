source("packages.R")
source("functions.R")
library(targets)


# undebug(plot_seu_marker_heatmap_by_scna)

tar_load(c("debranched_seus_1q", "debranched_seus_2p", "debranched_seus_6p", "debranched_seus_16q", "overall_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons"))

clone_dist_plots_2p <- map(unlist(debranched_seus_2p), ~plot_seu_marker_heatmap_by_scna_ara(.x, overall_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", return_plots = "clone_distribution_plot", h = 3, w = 5))

clone_dist_plots_6p <- map(unlist(debranched_seus_6p), ~plot_seu_marker_heatmap_by_scna_ara(.x, overall_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p", return_plots = "clone_distribution_plot", h = 3, w = 5))

collages <- map(unlist(debranched_seus_2p), ~plot_seu_marker_heatmap_by_scna_ara(.x, overall_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p"))


collages <- imap(collages, ~(.x + labs(title = .y)))

test0 <- wrap_plots(collages) + 
	plot_layout(ncol = 2) + 
	plot_annotation(title = "1q+")

down_pdf(test0)

down_pdf(test0[[1]], h = 12, w = 8)

down_pdf(test0[[2]], h = 12, w = 5)

down_pdf(test0[[6]], h = 3.5, w = 5)

