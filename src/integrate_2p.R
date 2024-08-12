library(patchwork)
library(targets)
source("packages.R")
source("functions.R")


# plot_seu_marker_heatmap_by_scna ------------------------------

tar_load(c("cluster_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons"))

seu_path <- "output/seurat/integrated_2p/integrated_seu_2p_complete.rds"
seu <- readRDS(seu_path)

seus <- SplitObject(seu, split.by = "batch")

seu_duo <-
	seus[c("SRR13884248", "SRR17960484")] |> 
	integration_workflow() |> 
	find_all_markers(seurat_assay = "integrated")

seu_duo$scna <-
	factor(ifelse(str_detect(seu_duo$GT_opt, "2"), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))

saveRDS(seu_duo, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")

seu_duo <- readRDS("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")

test1 <- make_clone_distribution_figure_debug(seu_duo, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", cluster_order = cluster_orders, group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p") |> 
	browseURL()

seu_duo@meta.data <- seu_duo@meta.data |> 
	tibble::rownames_to_column("cell") |> 
	dplyr::mutate(idents = 
									dplyr::case_when(`integrated_snn_res.0.4` %in% c(2,1,4) ~ "2p_enrich",
																	 `integrated_snn_res.0.4` %in% c(9) ~ "2p_depleted")) |> 
	tibble::column_to_rownames("cell") |> 
	identity()

mydiffex <- FindMarkers(seu_duo, group.by = "idents", ident.1 = "2p_enrich", ident.2 = "2p_depleted") 

enrichment_analysis(mydiffex) |> 
	plot_enrichment()



# debug(make_clone_distribution_figure)

# seu_duo0 <- seu_duo[,seu_duo$percent.mt < 3] |> 
# 	seurat_cluster(resolution = 0.4) |> 
# 	find_all_markers(metavar = "integrated_snn_res.0.4")
# 
# seu_duo2 <- drop_mt_cluster(seu_duo, group.by = "integrated_snn_res.0.4") |> 
# 	# seurat_cluster(resolution = 0.4) |> 
# 	find_all_markers(metavar = "integrated_snn_res.0.4") |>
# 	identity()
# 
# seu_duo3 <- drop_mt_cluster(seu_duo, group.by = "integrated_snn_res.0.4") |> 
# 	seurat_cluster(resolution = 0.4) |>
# 	find_all_markers(metavar = "integrated_snn_res.0.4") |>
# 	identity()
# 
# debug(make_clone_distribution_figure_debug)
# 
# test1 <- make_clone_distribution_figure_debug(seu_duo, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p")
# 
# test10 <- make_clone_distribution_figure_debug(seu_duo0, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p")
# 
# test20 <- make_clone_distribution_figure_debug(seu_duo2, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p") |> 
# 	browseURL()
# 
# test30 <- make_clone_distribution_figure_debug(seu_duo3, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p") |> 
# 	browseURL()

# test1 <- make_clone_distribution_figure("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", cluster_order = cluster_orders, group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p")


saveRDS(seu, seu_path)

debug(make_clone_distribution_figure)

test2 <- make_clone_distribution_figure(seu_path, group.bys = "clusters", scna_of_interest = "2p")

test0 <- make_clone_distribution_figure(seu_path, group.bys = "integrated_snn_res.0.4")

debug(make_clone_distribution_figure)




debug(assign_phase_clusters)

seu |> 
assign_phase_clusters(file_id = fs::path_file(seu_path), cluster_orders = cluster_orders) |> 
SplitObject(split.by = "batch") |> 
imap(~(saveRDS(.x, glue("{fs::path_dir(seu_path)}/{.y}_filtered_seu.rds")))) |>
identity()


test0 <- dir_ls(fs::path_dir(seu_path)) |> 
	set_names() |> 
	map(plot_clone_cc_plots, var_y = "clusters")


test0 <- test0[c(4,1:3)]

wrap_plots(test0) + 
	plot_layout(ncol = 1, guides= "collect") + 
	plot_annotation(title = "asdf")

plot_path <- ggsave("results/clone_cc_plots_16q.pdf", w = 6, h = 12)

browseURL(plot_path)

