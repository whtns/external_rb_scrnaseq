library(patchwork)
library(targets)
source("packages.R")
source("functions.R")


# plot_seu_marker_heatmap_by_scna ------------------------------

tar_load(c("cluster_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons"))

seu <- readRDS("output/seurat/integrated_16q/integrated_seu_16q_complete.rds")

# debug(assign_phase_clusters)

test0 <-
	seu |> 
	assign_phase_clusters(sample_id = "integrated_seu_16q_complete.rds", tumor_id = "", cluster_orders = cluster_orders) |> 
	SplitObject(split.by = "batch") |> 
	imap(~(saveRDS(.x, glue("output/seurat/integrated_16q/{.y}_filtered_seu.rds")))) |>
	identity()
	
	
	test0 <- dir_ls("output/seurat/integrated_16q/") |> 
		set_names() |> 
		map(plot_clone_cc_plots, var_y = "clusters")
	
	test0 <- test0[c(4,1:3)]
	
	wrap_plots(test0) + 
		plot_layout(ncol = 1, guides= "collect") + 
		plot_annotation(title = "asdf")
	
	plot_path <- ggsave("results/clone_cc_plots_16q.pdf", w = 6, h = 12)
	
	browseURL(plot_path)
	
	