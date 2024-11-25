library(targets)
source("packages.R")
source("functions.R")
library(plotgardener)

# plot_seu_marker_heatmap_by_scna ------------------------------

tar_load(c("debranched_seus_6p", "cluster_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons"))

debug(plot_seu_marker_heatmap_by_scna)

plot_seu_marker_heatmap_by_scna(debranched_seus_6p[[1]], cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p")

# fig_s03c------------------------------

# debug(plot_fig_s03c)
# debug(make_rb_scna_ideograms)
# debug(make_annoHighlight_from_consensus)

nb_paths <- dir_ls("output/numbat_sridhar/", regexp = ".*SRR[0-9]*_numbat.rds", recurse = TRUE) |> 
	sort()

make_rb_scna_ideograms(nb_paths[[3]]) |>
	# browseURL() |> 
	identity()

test0 <- plot_fig_s03c()

# find_diffex_clusters_between_corresponding_states ------------------------------

# debug(find_diffex_clusters_between_corresponding_states)

tar_load(c("states_dictionary_2p", "large_clone_comparisons", "numbat_rds_files"))

# test0 <- find_diffex_clusters_between_corresponding_states("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", states_dictionary_2p, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")

test0 <- readRDS("~/tmp.rds")

# plot_corresponding_clusters_diffex_volcanos 
# debug(plot_corresponding_clusters_diffex_volcanos)
# debug(make_volcano_plots)

test1 <- plot_corresponding_clusters_diffex_volcanos(test0, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds") |> 
	browseURL()

test2 <- plot_corresponding_enrichment(test0, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")


# effect of regression ------------------------------

# debug(plot_effect_of_regression)

plot_effect_of_regression("output/seurat/SRR14800534_filtered_seu.rds",
													"output/seurat/SRR14800534_regressed_seu.rds",
													w = 14, h = 14, filter_dropped_cluster = 7, regress_dropped_cluster = 7) |> 
	browseURL()


plot_effect_of_regression(final_seus, regressed_seus, w = 18, h = 12)






# plot_effect_of_filtering_regression ------------------------------

# tar_load("table_s03")

# debug(plot_fig_02_01_1)

table_s03 <- read_csv("results/table_s03.csv")

plot_fig_02_01_1("output/seurat/SRR13884249_unfiltered_seu.rds", "output/seurat/SRR13884249_filtered_seu.rds", removed_clusters = table_s03) |> 
	browseURL()

unfiltered_seu <- readRDS("output/seurat/SRR13884249_unfiltered_seu.rds")
filtered_seu <- readRDS("output/seurat/SRR13884249_filtered_seu.rds")
# str_extract(unfiltered_seu_path, "SRR[0-9]*")





# make oncoprint plots ------------------------------

tar_load(c("oncoprint_input_by_scna", "debranched_clone_trees", "oncoprint_settings"))

# debug(make_oncoprint_plots)

# debug(plot_recurrence)

make_oncoprint_plots(oncoprint_input_by_scna, debranched_clone_trees, oncoprint_settings, label = "_by_clone")


# make_clone_distribution_figure ------------------------------

# debug(make_clone_distribution_figure)

tar_load(c("cluster_orders"))

test2 <- make_clone_distribution_figure("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", cluster_order = cluster_orders, height = 12, width = 16) |> 
	browseURL()

test1 <- make_clone_distribution_figure("output/seurat/integrated_16q/integrated_seu_16q_complete.rds", cluster_order = cluster_orders, height = 12, width = 16) |> 
	browseURL()



# test0 <- make_clone_distribution_figure("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", cluster_order = cluster_orders)

# fig_02 ------------------------------

# debug(plot_fig_03)
# debug(plot_clone_cc_plots)
# debug(make_clone_distribution_figure)
# debug(plot_distribution_of_clones_across_clusters)

tar_load("cluster_orders")

# fig_02 <- plot_fig_02(cluster_orders)

fig_03 <- plot_fig_03(cluster_orders)

# make_clone_distribution_figure ------------------------------

# debug(make_clone_distribution_figure)

tar_load("cluster_orders")

make_clone_distribution_figure("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", cluster_order = cluster_orders)

make_clone_distribution_figure("output/seurat/integrated_1q/integrated_seu_1q_trio.rds", cluster_order = cluster_orders)

debug(find_diffex_clones)

tar_load(c("integrated_seus_1q", "numbat_rds_files", "large_clone_comparisons"))

diffex_clones_1q <- function(seu_path){
	seu <- readRDS(seu_path)
	
	diffex <- FindMarkers(seu, group.by = "scna", ident.1 = "w_scna", ident.2 = "wo_scna")
	
	excluded_genes <- 
		c(
			unlist(Seurat::cc.genes.updated.2019),
			str_subset(rownames(diffex), "^RP")
		)
	
	diffex_no_cc_no_rp <- diffex[!rownames(diffex) %in% excluded_genes,]
	
}

find_diffex_clones("output/seurat/integrated_1q/SRR13884249_integrated_1q_filtered_seu.rds", numbat_rds_files, large_clone_comparisons, location = "in_segment")




# make_clustree_for_speckle_set ------------------------------

# # debug(make_clustrees_for_sample)
# debug(make_clustree_for_clone_comparison)
# # # debug(plot_clustree_per_comparison)
# # debug(color_clustree_by_clone)
# debug(chi_sq_daughter_clusters)

# debug(plot_clustree_per_comparison)
# debug(color_clustree_by_clone)

test1 <- make_clustrees_for_sample("output/seurat/integrated_16q/integrated_seu_16q_complete.rds", mylabel = "16q_integrated", assay = "integrated")

test0 <- make_clustrees_for_sample("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", mylabel = "1q_integrated", assay = "integrated")

test2 <- make_clustrees_for_sample("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", mylabel = "2p_integrated", assay = "integrated")


# make_numbat_heatmaps ------------------------------

tar_load(c("original_seus", "numbat_rds_files"))

# debug(make_numbat_heatmaps)
debug(plot_numbat)
# debug(plot_phylo_heatmap)
debug(plot_variability_at_SCNA)

make_numbat_heatmaps(original_seus[[1]], numbat_rds_files, p_min = 0.5, line_width = 0.1, extension = "_filtered") |> 
	map(browseURL)


# plot_fig_09_10 ------------------------------

tar_load(c("corresponding_seus_2p", "corresponding_seus", "corresponding_clusters_diffex", "corresponding_clusters_enrichments"))

debug(plot_fig_09_10)
debug(compare_corresponding_enrichments)

plot_fig_09_10(corresponding_seus_2p, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments, recurrence_threshold = 2, plot_path = "results/fig_09.pdf", widths = rep(4, 3), heights = c(8,4,8), common_seus = c("SRR13884248_filtered_seu_2p.rds", "SRR17960484_filtered_seu_2p.rds"))


debug(plot_fig_s20)

test3 <- plot_fig_s20()

# make_table_s03 ------------------------------

make_table_s03()

# debug(plot_fig_s20)

fig_s20 <- plot_fig_s20()

# debug(plot_fig_s09)

fig_s09 <- plot_fig_s09()

# make_table_s07_s09 ------------------------------

# debug(make_table_s07_s09)

test2 <- make_table_s10(seu_path = "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", table_path = "results/table_s10.csv")

test0 <- make_table_s07(seu_path = "output/seurat/integrated_1q/integrated_seu_1q_complete.rds", table_path = "results/table_s07.csv")

test1 <- make_table_s09(seu_path = "output/seurat/integrated_16q/integrated_seu_16q_complete.rds", table_path = "results/table_s09.csv")

# plot fig_s06 ------------------------------

# debug(plot_fig_s06)

fig_s06 <- plot_fig_s06()

# plot_fig_07c ------------------------------

# debug(plot_fig_07c)
plot_fig_07c()

# find_diffex_bw_clones_for_each_cluster_integrated ------------------------------

tar_load(c("debranched_seus_1q", "integrated_seus_1q", "integrated_seus_16q", "numbat_rds_files",  "large_clone_comparisons", "debranched_cluster_orders"))

# debug(find_diffex_bw_clones_for_each_cluster_integrated)
# debug(make_clone_comparison)
# debug(clone_diff_per_cluster)

test0 <- find_diffex_bw_clones_for_each_cluster_integrated("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", numbat_rds_files, large_clone_comparisons, location = "all", scna_of_interest = "1q+")



# find_diffex_bw_clones_for_each_cluster ------------------------------

tar_load(c("debranched_seus_1q", "integrated_seus_1q", "integrated_seus_16q", "numbat_rds_files",  "large_clone_comparisons", "debranched_cluster_orders"))

debug(find_diffex_bw_clones_for_each_cluster)
debug(make_clone_comparison)
debug(clone_diff_per_cluster)

test0 <- find_diffex_bw_clones_for_each_cluster(integrated_seus_16q[[1]], numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-")

# plot_fig_07_08 ------------------------------

tar_load('fig_07b_input')

debug(plot_fig_07_08)

plot_fig_07_08(fig_07b_input, plot_path = "results/fig_07.pdf", height = 5, width = 4)




# plot_fig_07_08 ------------------------------

tar_load(c("fig_07a_input", "fig_08_input"))

# debug(plot_fig_07_08)

test0 <- map(c(0.99, 0.5, 0.1, 0.05), ~plot_fig_07_08(fig_07a_input, plot_title = .x, p_adj_threshold = .x, plot_path = glue("results/fig_07_{.x}.pdf"), w = 8, h = 8)) |> 
	map(browseURL)

test0 <- map(c(0.99, 0.5, 0.1, 0.05), ~plot_fig_07_08(fig_08_input, plot_title = .x, p_adj_threshold = .x, plot_path = glue("results/fig_08_{.x}.pdf"), w = 8, h = 8)) |> 
	map(browseURL)



# plot_study_metadata ------------------------------

tar_load("study_cell_stats")

# debug(plot_study_metadata)
# debug(plot_mt_v_nUMI)

test0 <- plot_study_metadata(study_cell_stats)

# find_diffex_bw_clones_for_each_cluster ------------------------------

# debug(find_diffex_bw_clones_for_each_cluster)
# debug(clone_diff_per_cluster)

tar_load(c("integrated_seus_1q", "numbat_rds_files",  "large_clone_comparisons", "debranched_cluster_orders"))

debug(find_diffex_bw_clones_for_each_cluster)
debug(make_clone_comparison)
debug(clone_diff_per_cluster)

find_diffex_bw_clones_for_each_cluster(integrated_seus_1q[[1]], numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all")




# plot_seu_marker_heatmap_all_resolutions ------------------------------

debug(plot_seu_marker_heatmap_all_resolutions)
debug(plot_seu_marker_heatmap)
# debug(plot_clone_tree)

tar_load(c("debranched_seus_6p", "debranched_cluster_orders", "large_clone_simplifications", "numbat_rds_files"))

test0 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[10]], numbat_rds_files, large_clone_simplifications)

test1 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[11]], numbat_rds_files, large_clone_simplifications, debranched_cluster_orders)

test2 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[12]], numbat_rds_files, large_clone_simplifications, debranched_cluster_orders)

# filter_cluster_save_seu ------------------------------

tar_load(c("numbat_rds_files", "cluster_dictionary", "large_clone_simplifications", "cells_to_remove"))

debug(filter_cluster_save_seu)

filter_cluster_save_seu("output/numbat_sridhar/SRR14800534_numbat.rds", cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove, extension = "", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv")


# find_diffex_clusters_between_corresponding_states ------------------------------

tar_load(c("corresponding_seus", "corresponding_states_dictionary", "large_clone_comparisons", "numbat_rds_files"))

# debug(find_diffex_clusters_between_corresponding_states)

find_diffex_clusters_between_corresponding_states(corresponding_seus[[8]], corresponding_states_dictionary[8], large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")


# plot_fig_09_10 ------------------------------

debug(plot_fig_09_10)
# debug(dotplot_recurrent_genes)
# debug(compare_corresponding_enrichments)

tar_load(c("corresponding_seus_2p", "corresponding_seus_6p", "corresponding_seus", "corresponding_clusters_diffex", "corresponding_clusters_enrichments"))

fig_10 <- plot_fig_09_10(corresponding_seus_6p, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments, recurrence_threshold = 2, plot_path = "results/fig_10.pdf", widths = rep(4,3), heights = c(12, 4, 12), common_seus = c("SRR13884248_filtered_seu_6p.rds", "SRR17960484_filtered_seu_6p.rds"))


# plot_fig_04_05 ------------------------------

tar_load(c("corresponding_seus_2p", "corresponding_seus", "corresponding_clusters_diffex", "corresponding_clusters_enrichments", "cluster_orders", "integrated_seus_2p", "integrated_seus_6p"))

# debug(plot_fig_04_05)

fig_04 <- plot_fig_04_05(c("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", "output/seurat/SRR13884247_branch_6_filtered_seu.rds"), corresponding_clusters_enrichments[[6]], integrated_seu_paths = integrated_seus_2p, cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", width = 18, height = 10)


fig_05 <- plot_fig_04_05("output/seurat/integrated_6p/integrated_seu_6p_complete.rds", corresponding_clusters_enrichments[[7]], integrated_seu_paths = integrated_seus_6p, plot_path = "results/fig_05.pdf", cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p", width = 18, height = 10)

# compare_corresponding_enrichments ------------------------------

debug(compare_corresponding_enrichments)

tar_load(c("corresponding_states_dictionary", "corresponding_seus", "large_clone_comparisons", "numbat_rds_files", "corresponding_clusters_diffex"))

test2 <- map2(corresponding_clusters_diffex, corresponding_seus, compare_corresponding_enrichments)



tar_load(c("cluster_orders", "debranched_seus_2p", "debranched_seus_6p", "numbat_rds_files", "large_clone_comparisons"))



# find_diffex_clones_between_corresponding_states ------------------------------

tar_load(c("corresponding_states_dictionary", "corresponding_seus", "large_clone_comparisons", "numbat_rds_files", "corresponding_clusters_diffex"))

# debug(find_diffex_clusters_between_corresponding_states)
# debug(make_cluster_comparison)

test0 <- find_diffex_clusters_between_corresponding_states(corresponding_seus[[6]], corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")

test5 <- find_diffex_clusters_between_corresponding_states(corresponding_seus[[5]], corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")





fig_03 <- plot_fig_03(cluster_orders)

fig_04 <- plot_fig_04(c("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", "output/seurat/SRR13884247_branch_6_filtered_seu.rds"), cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", width = 18, height = 10) |> 
	map(browseURL)

# debug(plot_fig_04_panels)

plot_fig_04(c("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", "output/seurat/SRR13884247_branch_6_filtered_seu.rds"), cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", width = 18, height = 18) |> 
	map(browseURL)

# plot_fig_04 ------------------------------

tar_load(c("cluster_orders", "debranched_seus_2p", "debranched_seus_6p", "numbat_rds_files", "large_clone_comparisons"))

# debug(plot_fig_04)
# debug(drop_mt_cluster)
# debug(seu_gene_heatmap)
# debug(seu_complex_heatmap2)

# test0 <- plot_fig_04(debranched_seus_2p[[3]], cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, nb_paths = numbat_rds_files, scna_of_interest = "2p", width = 16, height = 12) |> 
# 	browseURL()

plot_fig_04(seu_paths, cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, nb_paths = numbat_rds_files, scna_of_interest = "2p", width = 16, height = 12)

test1 <- plot_fig_04("output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds", cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, nb_paths = numbat_rds_files, scna_of_interest = "2p", width = 16, height = 12)

test2p <- map(debranched_seus_2p, plot_fig_04 , cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, nb_paths = numbat_rds_files, scna_of_interest = "2p", width = 16, height = 12)

test6p <- map(debranched_seus_6p, plot_fig_05 , cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, nb_paths = numbat_rds_files, scna_of_interest = "6p", width = 16, height = 12)


# plot_seu_marker_heatmap_by_scna ------------------------------

tar_load(c("cluster_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons"))

seu_path <- "output/seurat/integrated_2p/integrated_seu_2p_complete.rds"
seu <- readRDS(seu_path)

debug(make_clone_distribution_figure)

test1 <- make_clone_distribution_figure(seu_path, cluster_order = cluster_orders, group.bys = "integrated_snn_res.0.4", scna_of_interest = "2p")

# plot_corresponding_enrichment ------------------------------

tar_load(c("corresponding_clusters_diffex", "corresponding_seus"))

debug(plot_corresponding_enrichment)

test0 <- plot_corresponding_enrichment(corresponding_clusters_diffex[[1]], unlist(corresponding_seus)[[1]])

# plot_corresponding_clusters_diffex ------------------------------

tar_load(c("corresponding_clusters_diffex", "corresponding_seus", "corresponding_states_dictionary", "large_clone_comparisons", "numbat_rds_files"))

debug(plot_corresponding_enrichment)

test1 <- plot_corresponding_enrichment(corresponding_clusters_diffex[[1]], corresponding_seus[[1]], corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")

test0 <- map2(corresponding_clusters_diffex, corresponding_seus, plot_corresponding_clusters_diffex_volcanos)

test1 <- map2(corresponding_clusters_diffex, corresponding_seus, plot_corresponding_enrichment, corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")

# plot_clone_cc_plots ------------------------------

tar_load(c("debranched_seus_1q", "debranched_seus_2p", "debranched_seus_6p", "debranched_seus_16q", "cluster_orders", "large_clone_comparisons"))

integrated_1q_seus <- dir_ls("output/seurat/integrated_1q/", glob = "*.rds")

names(integrated_1q_seus) <- fs::path_ext_remove(fs::path_file(integrated_1q_seus))



clone_cc_plots_by_scna <- function(seu_list, scna_of_interest = "1q", ...) {
	test0 <- map(seu_list, plot_clone_cc_plots, scna_of_interest = scna_of_interest, ...)
	print(test0)
	
	myplot <- wrap_plots(test0) + 
		plot_layout(ncol = 1)
	
	plot_path <- ggsave(glue("results/{scna_of_interest}_clone_distribution.pdf"), myplot, h = length(test0)*16/5, w = length(test0)*8/5)
	
	return(plot_path)
}

test0 <- clone_cc_plots_by_scna(debranched_seus_16q, scna_of_interest = "16q", large_clone_comparisons = large_clone_comparisons)

browseURL(test0)




# make_expression_heatmap_comparison ------------------------------

tar_load(c("large_numbat_pdfs", "filtered_numbat_heatmaps"))

debug(make_expression_heatmap_comparison)

make_expression_heatmap_comparison(large_numbat_pdfs, filtered_numbat_heatmaps)


# plot_seu_marker_heatmap_by_scna ------------------------------

tar_load(c("cluster_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons", "debranched_seus_2p", "debranched_seus_6p"))


# make_cc_plot(integrated_1q_seu, "integrated_snn_res.0.4")
# 
# make_cc_plot(integrated_1q_seu, "integrated_snn_res.0.6")

# debug(make_clone_distribution_figure)

# debug(find_cluster_pairwise_distance)

map(debranched_seus_2p, make_clone_distribution_figure, scna_of_interest = "2p\\+", width =20, height = 10, group.bys = c("clusters")) |>
	map(browseURL)

map(debranched_seus_6p, make_clone_distribution_figure, scna_of_interest = "6p\\+", width =20, height = 10, group.bys = c("clusters")) |>
	map(browseURL)

# 2p ------------------------------

scna_of_interest = "2p\\+"
for(seu_path in debranched_seus_2p){
	new_colnames <- glue("new_SCT_snn_res.{seq(0.2, 2.0, by = 0.2)}") |> 
		set_names()
	
	old_colnames <- glue("old_SCT_snn_res.{seq(0.2, 2.0, by = 0.2)}") |> 
		set_names()
	
	# debug(make_clone_distribution_figure)
	
	make_clone_distribution_figure(seu_path, scna_of_interest = scna_of_interest, width =18, height = 10, group.bys = c(rbind(new_colnames, old_colnames)))
}

# 6p ------------------------------

scna_of_interest = "6p\\+"
for(seu_path in debranched_seus_6p){
	new_colnames <- glue("new_SCT_snn_res.{seq(0.2, 2.0, by = 0.2)}") |> 
		set_names()
	
	old_colnames <- glue("old_SCT_snn_res.{seq(0.2, 2.0, by = 0.2)}") |> 
		set_names()
	
	# debug(make_clone_distribution_figure)
	
	make_clone_distribution_figure(seu_path, scna_of_interest = scna_of_interest, width =18, height = 10, group.bys = c(rbind(new_colnames, old_colnames)))
}


# end set ------------------------------

map(debranched_seus_2p, make_clone_distribution_figure, scna_of_interest = "2p\\+", width =18, height = 10, group.bys = c("clusters"))

map(debranched_seus_6p, make_clone_distribution_figure, scna_of_interest = "6p\\+", width =18, height = 10, group.bys = c("clusters"))


make_clone_distribution_figure("output/seurat/integrated_6p/integrated_seu_6p_complete.rds", cluster_order = cluster_orders, scna_of_interest = "6p", width =14, height = 10) |> 
	browseURL()

make_clone_distribution_figure("output/seurat/integrated_16q/integrated_seu_16q_complete.rds", cluster_order = cluster_orders, scna_of_interest = "16q", width =14, height = 10) |> 
	browseURL()

make_clone_distribution_figure("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", cluster_order = cluster_orders, 
															 scna_of_interest = "1q", width =18, height = 16) |> 
	browseURL()

debug(plot_seu_marker_heatmap_by_scna)

plot_seu_marker_heatmap_by_scna("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", cluster_orders, numbat_rds_files, large_clone_simplifications = large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q") |> 
	browseURL()

plot_seu_marker_heatmap_by_scna("output/seurat/SRR13884249_filtered_seu_1q_integrated.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q") |> 
	browseURL()

tar_load(c("numbat_rds_files", "large_clone_simplifications", "cluster_orders", "rb_scna_samples", "large_clone_comparisons"))

plot_seu_marker_heatmap_by_scna("output/seurat/SRR14800534_filtered_seu.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q") |> 
	browseURL()

plot_seu_marker_heatmap_by_scna("output/seurat/SRR14800535_filtered_seu_1q_integrated.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q") |> 
	browseURL()

plot_seu_marker_heatmap_by_scna("output/seurat/SRR14800536_filtered_seu_1q_integrated.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q") |> 
	browseURL()

plot_seu_marker_heatmap_by_scna("output/seurat/SRR14800534_filtered_seu_16q_integrated.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q")

plot_seu_marker_heatmap_by_scna("output/seurat/SRR14800535_filtered_seu_16q_integrated.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q")

plot_seu_marker_heatmap_by_scna("output/seurat/SRR14800536_filtered_seu_16q_integrated.rds", cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q")





# plot_fig_03_collage ------------------------------

tar_load(c("debranched_seus_2p", "cluster_orders", "numbat_rds_files", "large_clone_simplifications", "rb_scna_samples", "large_clone_comparisons"))

plot_fig_04(debranched_seus_2p[[1]], cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p")


browseURL(plot_path)

# plot_corresponding_clusters_diffex_heatmaps ------------------------------

tar_load(c("corresponding_states_dictionary", "corresponding_seus", "large_clone_comparisons", "numbat_rds_files", "corresponding_clusters_diffex"))

# debug(plot_corresponding_clusters_diffex_heatmaps)

test0 <- 
	plot_corresponding_clusters_diffex_heatmaps(corresponding_clusters_diffex[[1]], corresponding_seus[[1]], corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")

# plot_corresponding_clusters_diffex_volcanos ------------------------------

tar_load(c("corresponding_states_dictionary", "corresponding_seus", "large_clone_comparisons", "numbat_rds_files", "corresponding_clusters_diffex"))

# debug(plot_corresponding_clusters_diffex_volcanos)

test0 <- map2(corresponding_clusters_diffex, corresponding_seus, plot_corresponding_clusters_diffex_volcanos)




# plot_clone_pearls ------------------------------

tar_load(c("debranched_seus_1q", "debranched_seus_2p", "debranched_seus_6p", "debranched_seus_16q", "cluster_orders", "large_clone_comparisons"))

# debug(plot_clone_pearls)
# debug(plot_distribution_of_clones_pearls)
debug(make_pearls_plot)

test0 <- plot_clone_pearls(debranched_seus_2p[[2]], var_y = "clusters")

# enrich_by_cluster ------------------------------

tar_load(c("debranched_seus_2p"))

test0 <- enrich_by_cluster(debranched_seus_2p[[1]])



tar_load(c("debranched_seus_1q", "debranched_seus_2p", "debranched_seus_6p", "debranched_seus_16q", "cluster_orders", "large_clone_comparisons", "cluster_orders"))

# debug(calculate_clone_distribution)

names(debranched_seus_1q) <- fs::path_file(debranched_seus_1q)

# debug(make_pairwise_plots)

clone_dist_res_1q <- map(debranched_seus_1q, calculate_clone_distribution, cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q", pairwise = TRUE)

clone_dist_res_2p <- map(debranched_seus_2p, calculate_clone_distribution, cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", pairwise = TRUE)

clone_dist_res_6p <- map(debranched_seus_6p, calculate_clone_distribution, cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p", pairwise = TRUE)

clone_dist_res_16q <- map(debranched_seus_16q, calculate_clone_distribution, cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q", pairwise = TRUE)

pull_chosen_cluster_resolution <- function(clone_dist_res, resolution_dictionary, scale = "absolute"){
	
	table_excel_paths <- map(clone_dist_res, "table")
	
	selected_sheets <- resolution_dictionary[names(table_excel_paths)]
	
	phase_level_df <- map2(table_excel_paths, selected_sheets, ~{readxl::read_excel(.x, sheet = .y)}) |> 
		map(set_names, c("clusters", "wo_scna", "w_scna", "sample_id", "up", "down")) |> 
		dplyr::bind_rows(.id = "file_id") |>
		tibble::rownames_to_column("paired") |> 
		tidyr::pivot_longer(cols = !any_of(c("clusters", "sample_id", "up", "down", "file_id", "paired")), names_to = "scna_status", values_to = "percent") |>
		identity()
	
	if(scale == "relative"){
		phase_level_df <- phase_level_df
	}
	
	myplot <- 
		phase_level_df |> 
		dplyr::mutate(clusters = factor(clusters, levels = c("all", "g1", "g1_s", "s", "s_g2", "g2_m", "pm", "hypoxia", "hsp", "other", "s_star"))) |> 
		ggplot(aes(x = scna_status, y = percent, fill = scna_status)) +
		geom_line(aes(group = paired), size=0.5) +
		geom_boxplot()+
		geom_point(position = position_dodge(width = 0.75)) +
		facet_wrap(~clusters, scales = "free_x", nrow = 1) +
		theme_minimal() + 
		theme(
			axis.text.x = element_blank()
		)
		NULL
	
	return(myplot)
	
}

pdf("results/groups_clone_proportions.pdf", w= 8, h = 4)
pull_chosen_cluster_resolution(clone_dist_res_1q, resolution_dictionary, scale = "relative") +
	labs(title = "1q+ clone proportions")

pull_chosen_cluster_resolution(clone_dist_res_2p, resolution_dictionary) + 
	labs(title = "2p+ clone proportions")

pull_chosen_cluster_resolution(clone_dist_res_6p, resolution_dictionary) + 
	labs(title = "6p+ clone proportions")

pull_chosen_cluster_resolution(clone_dist_res_16q, resolution_dictionary) + 
	labs(title = "16q- clone proportions")
dev.off()

browseURL("results/groups_clone_proportions.pdf")

# assign_designated_phase_clusters ------------------------------

tar_load(c("scna_seus", "cluster_orders", "resolution_dictionary"))

assign_designated_phase_clusters(scna_seus[[2]], cluster_orders, resolution_dictionary)

# plot_seu_gene_heatmap ------------------------------

tar_load(c("large_clone_comparisons", "debranched_seus_1q"))

# debug(plot_seu_gene_heatmap)
# debug(seu_gene_heatmap)
# debug(seu_complex_heatmap2)

test0 <- plot_seu_gene_heatmap(debranched_seus_1q[[4]] , large_clone_comparisons, scna_of_interest = "1q", w = 8, h = 12)

test0 <- map(debranched_seus_1q, plot_seu_gene_heatmap , large_clone_comparisons, scna_of_interest = "1q", w = 8, h = 12)

test0 <- map(debranched_seus_2p, plot_seu_gene_heatmap , large_clone_comparisons, scna_of_interest = "2p", w = 8, h = 12)

test0 <- map(debranched_seus_6p, plot_seu_gene_heatmap , large_clone_comparisons, scna_of_interest = "6p", w = 8, h = 12)

test0 <- map(debranched_seus_16q, plot_seu_gene_heatmap , large_clone_comparisons, scna_of_interest = "16q", w = 8, h = 12)


plot_putative_marker_across_samples(interesting_genes[["S2"]], debranched_seus_16q, plot_type = VlnPlot, group_by = "clusters", extension = "s2")

plot_putative_marker_across_samples(interesting_genes, scna_seus, plot_type = VlnPlot, group_by = "scna", cluster_dictionary)

debug(assign_designated_phase_clusters)
debug(assign_phase_cluster_at_resolution)

assign_designated_phase_clusters(debranched_seus_2p[[6]], cluster_orders, resolution_dictionary)

seu_path <- debranched_seus_2p[[6]]

# marker_gene_vlnplots_by_cluster ------------------------------

tar_load(c("interesting_genes", "debranched_seus_1q"))

# debug(plot_markers_in_sample)

plot_putative_marker_across_samples(interesting_genes, unlist(debranched_seus_1q), plot_type = VlnPlot, group_by = "clusters", cluster_dictionary)



# plot_seu_marker_heatmap ------------------------------

tar_load(c("debranched_seus", "debranched_cluster_orders", "large_clone_simplifications", "numbat_rds_files"))

debug(plot_seu_marker_heatmap)
# debug(plot_clone_tree)
debug(seu_complex_heatmap)

# plot_seu_marker_heatmap(filtered_seus[[1]], cluster_orders[[1]])

# plot_seu_marker_heatmap(seu_path = filtered_seus[[1]], nb_path = numbat_rds_files[[1]], clone_simplifications = large_clone_simplifications)

# test0 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[13]], numbat_rds_files, large_clone_simplifications)

test0 <- plot_seu_marker_heatmap(debranched_seus[[12]], debranched_cluster_orders, numbat_rds_files, large_clone_simplifications)

browseURL(test0)

plot_seu_marker(debranched_seus, debranched_cluster_orders, numbat_rds_files, large_clone_simplifications)



# plot_study_metadata ------------------------------

tar_load("study_cell_stats")
# 
# debug(plot_study_metadata)
# debug(plot_study_cell_stats)

test0 <- plot_study_metadata(study_cell_stats)



# filter_oncoprint_diffex ------------------------------

tar_load(c("unfiltered_oncoprint_input_by_scna", "oncoprint_settings"))
# debug(filter_oncoprint_diffex)

# debug(filter_input_by_scna)

oncoprint_input_by_scna <- filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna, oncoprint_settings)

tar_load(c("oncoprint_settings", "debranched_clone_trees"))

# debug(make_oncoprint_plots)
# debug(plot_recurrence)

test0 <- make_oncoprint_plots(oncoprint_input_by_scna, debranched_clone_trees, oncoprint_settings, label = "_by_clone")

# make_oncoprint_plots(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans, oncoprint_input_by_scna$all, debranched_clone_trees, oncoprint_settings, label = "_by_clone", p_val_threshold = 0.1)




# calculate_clone_distribution ------------------------------

tar_load(c("debranched_seus", "debranched_cluster_orders"))

# debug(calculate_clone_distribution)
# debug(make_pairwise_plots)
# debug(plot_distribution_of_clones_across_clusters)

calculate_clone_distribution(debranched_seus[[10]], cluster_orders, pairwise = TRUE)


# plot_distribution_of_clones_across_clusters ------------------------------

tar_load("debranched_seus")

debug(plot_distribution_of_clones_across_clusters)

seu <- readRDS(debranched_seus[[1]])

plot_distribution_of_clones_across_clusters(seu, "test")


# make_oncoprint_plots for each cluster------------------------------

tar_load(c("oncoprint_input_by_scna_for_each_cluster", "debranched_clone_trees", "oncoprint_settings_by_cluster"))

debug(make_oncoprint_plots)
# debug(plot_recurrence)

make_oncoprint_plots(oncoprint_input_by_scna_for_each_cluster$cis, oncoprint_input_by_scna_for_each_cluster$trans, oncoprint_input_by_scna_for_each_cluster$all, debranched_clone_trees, oncoprint_settings_by_cluster, label = "_by_cluster", p_val_threshold = 1)


# make_oncoprint_diffex ------------------------------

tar_load(c("large_filter_expressions", "cluster_dictionary", "debranched_ids", "cis_diffex_clones_for_each_cluster", "trans_diffex_clones_for_each_cluster", "all_diffex_clones_for_each_cluster", "large_clone_comparisons", "rb_scna_samples"))

debug(make_oncoprint_diffex)

make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, all_diffex_clones_for_each_cluster, large_clone_comparisons, rb_scna_samples, by_cluster = TRUE, n_slice = 20)



# make_oncoprint_plots ------------------------------

tar_load(c("oncoprint_input_by_scna", "debranched_clone_trees"))

debug(make_oncoprint_plots)
# debug(plot_recurrence)

make_oncoprint_plots(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans, oncoprint_input_by_scna$all, debranched_clone_trees)

# compile_cis_trans_enrichment_recurrence ------------------------------

tar_load("oncoprint_enrich_clones_hallmark")

# debug(compile_cis_trans_enrichment_recurrence)
# debug(plot_enrichment_recurrence)

test1 <- 
	compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_hallmark,
																					cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_hallmark.pdf",
																					trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_hallmark.pdf",
																					cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_hallmark.xlsx",
																					trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_hallmark.xlsx"
	)

tar_load("oncoprint_enrich_clusters_hallmark")

test0 <- compile_cis_trans_enrichment_recurrence_by_cluster(oncoprint_enrich_clusters_hallmark,
																				cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_hallmark.pdf",
																				trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_hallmark.pdf",
																				cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_hallmark.xlsx",
																				trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_hallmark.xlsx",
																				by_cluster = TRUE
)




# debranch_seus------------------------------

tar_load(c("final_seus", "branch_dictionary", "debranched_cluster_orders", "numbat_rds_files", "large_clone_simplifications"))

debug(debranch_seus)
debug(split_seu_by_branch)
debug(prep_seu_branch)
debug(assign_phase_clusters)

debranch_seus(final_seus[3], 
							branch_dictionary, 
							cluster_orders = debranched_cluster_orders, 
							nb_paths = numbat_rds_files, 
							clone_simplifications = large_clone_simplifications)


# plot_seu_marker_heatmap_integrated ------------------------------

tar_load(integrated_seus)

debug(plot_seu_marker_heatmap_integrated)

test0 <- plot_seu_marker_heatmap_integrated(integrated_seus[[1]], group.by = "integrated_snn_res.0.2", assay = "integrated")
browseURL(test0)



# enrich_oncoprints_clusters ------------------------------

tar_load(c("cluster_dictionary", "interesting_samples", "cis_diffex_clones_for_each_cluster", "trans_diffex_clones_for_each_cluster", "large_clone_comparisons")
)

debug(enrich_oncoprints_clusters)


enrich_oncoprints_clusters(large_filter_expressions,
													 cluster_dictionary,
													 debranched_ids,
													 cis_diffex_clones_for_each_cluster,
													 trans_diffex_clones_for_each_cluster,
													 large_clone_comparisons,
													 gene_set = "hallmark"
)




# make_cluster_comparisons_by_phase_for_disctinct_clones ------------------------------

tar_load(c("overall_seus", "cluster_comparisons"))

debug(make_cluster_comparisons_by_phase_for_disctinct_clones)

debug(table_and_plot_enrichment)

make_cluster_comparisons_by_phase_for_disctinct_clones(cluster_comparisons[2], overall_seus)


# find_diffex_clones ------------------------------

debug(find_diffex_clones)
debug(make_clone_comparison)

tar_load(c("large_clone_comparisons", "numbat_rds_files", "debranched_seus"))

find_diffex_clones(debranched_seus[[22]], numbat_rds_files, large_clone_comparisons, location = "out_of_segment")


# enrich_oncoprints ------------------------------


tar_load(c(
	"large_filter_expressions",
	"cluster_dictionary",
	"debranched_ids",
	"cis_diffex_clones",
	"trans_diffex_clones",
	"all_diffex_clones",
	"large_clone_comparisons"
))

# debug(enrich_oncoprints_clusters)
debug(enrich_oncoprints)
# debug(enrichment_analysis)
# debug(compile_cis_trans_enrichment_recurrence)


enrich_oncoprints(large_filter_expressions,
									cluster_dictionary,
									debranched_ids,
									cis_diffex_clones,
									trans_diffex_clones,
									all_diffex_clones,
									large_clone_comparisons,
									gene_set = "hallmark"
)



# make_oncoprint_diffex ------------------------------

tar_load(c("large_filter_expressions", "cluster_dictionary", "debranched_ids", "cis_diffex_clones", "trans_diffex_clones", "all_diffex_clones", "large_clone_comparisons", "rb_scna_samples"))

debug(make_oncoprint_diffex)

make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, rb_scna_samples, n_slice = 20)



# assign_phase_clusters ------------------------------

debug(assign_phase_clusters)

tar_load(c("debranched_seus", "debranched_cluster_orders", "numbat_rds_files", "large_clone_simplifications"))

# map(debranched_seus, ~assign_phase_clusters(.x, debranched_cluster_orders, numbat_rds_files, large_clone_simplifications))

assign_phase_clusters(debranched_seus[[10]]	, debranched_cluster_orders, numbat_rds_files, large_clone_simplifications)



# pull_cluster_orders ------------------------------

# debug(pull_cluster_orders)

pull_cluster_orders("data/debranched_cluster_ids.csv")


# make_volcano_diffex_clones------------------------------

tar_load(c("cis_diffex_clones_for_each_cluster", "cis_diffex_clones"))

debug(make_volcano_diffex_clones)

make_volcano_diffex_clones(
	cis_diffex_clones_for_each_cluster,
	"results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
	cis_diffex_clones,
	"results/diffex_bw_clones_large_in_segment.pdf")



# make_numbat_heatmaps ------------------------------

tar_load("large_numbat_heatmap_files")

large_numbat_heatmap_files <- 
	large_numbat_heatmap_files %>% 
	set_names(str_extract(., "SRR[0-9]*"))


failed_numbat_samples <- c("SRR14800537")

compile_bad_scna_pdfs <- function(bad_scna_samples = c("SRR13884240", "SRR13884241", 
																												"SRR13884244", "SRR13884245", "SRR17960483")) {
  
  numbat_done_files <- dir_ls("output/numbat_sridhar/", glob = "*done.txt", recurse = TRUE) %>% 
  	set_names(str_extract(., "SRR[0-9]*"))
  
  
  panel_2_pngs <- numbat_done_files[names(numbat_done_files) %in% bad_scna_samples] %>% 
  	path_dir() %>% 
  	fs::path("panel_2.png") %>% 
  	set_names(str_extract(., "SRR[0-9]*"))
  
  panel_2_pdfs <- stringr::str_replace(panel_2_pngs, ".png", ".pdf")
  
  panel_2_images <- purrr::map(panel_2_pngs,image_read) %>%
  	purrr::imap(~image_annotate(.x, .y, size = 50))
  
  map2(panel_2_images, panel_2_pdfs, ~image_write(.x, format = "pdf", .y))
  
  bad_scna_colormap_compilation <- qpdf::pdf_combine(panel_2_pdfs, "results/bad_scnas.pdf")
}

test0 <- compile_bad_scna_pdfs()



make_numbat_heatmaps("output/seurat/SRR17960483_seu.rds", "output/numbat_sridhar/SRR17960483_numbat.rds", p_min = 0.5, line_width = 0.1, extension = "_unfiltered")


# convert_numbat_pngs ------------------------------

# debug(convert_numbat_pngs)

convert_numbat_pngs("output/numbat_sridhar/SRR13884242_numbat.rds")



# find_diffex_clones ------------------------------

tar_load(c("debranched_seus", "numbat_rds_files", "large_clone_comparisons"))

debug(find_diffex_clones)

find_diffex_clones(debranched_seus[[22]], numbat_rds_files, large_clone_comparisons, location = "in_segment")






# plot_clone_tree_from_path ------------------------------

debug(plot_clone_tree_from_path)

tar_load(c("debranched_seus", "numbat_rds_files", "large_clone_simplifications"))

plot_clone_tree_from_path(debranched_seus[[3]], numbat_rds_files, large_clone_simplifications, legend = FALSE, horizontal = FALSE)


# pull_cluster_orders ------------------------------
debug(pull_cluster_orders)

tar_load(debranched_cluster_file)

pull_cluster_orders(debranched_cluster_file,
										"data/debranched_cluster_ids.tsv")
 



# plot_clone_tree ------------------------------

tar_load("large_clone_simplifications")

debug(plot_clone_tree)

seu <- readRDS("output/seurat/SRR13884246_branch_6_filtered_seu.rds")

plot_clone_tree(seu, "SRR13884246", "output/numbat_sridhar/SRR13884246_numbat.rds", large_clone_simplifications)

# score_and_vlnplot_seu ------------------------------

tar_load(c("final_seus", "subtype_markers", "numbat_rds_files", "large_clone_simplifications"))

debug(score_and_vlnplot_seu)

score_and_vlnplot_seu(final_seus[3], numbat_rds_files[3], large_clone_simplifications[3], subtype_markers)



# prep_unfiltered_seu ------------------------------

tar_load(c("numbat_rds_files", "cluster_dictionary", "large_clone_simplifications", "large_filter_expressions", "cells_to_remove"))

# debug(prep_unfiltered_seu)

prep_unfiltered_seu(numbat_rds_files[15], cluster_dictionary, large_clone_simplifications, large_filter_expressions, extension = "_unfiltered")



# make_numbat_plot_files ------------------------------

tar_load(c("numbat_rds_files", "final_seus", "cluster_dictionary", "large_filter_expressions", "large_clone_simplifications"))

debug(make_numbat_plot_files)

make_numbat_plot_files(numbat_rds_files[[15]], final_seus[[15]], cluster_dictionary[[15]], large_filter_expressions[[15]], large_clone_simplifications, extension = "_filtered")

# collect_study_metadata ------------------------------

debug(collect_study_metadata)

collect_study_metadata()


# find_diffex_from_clustree ------------------------------

tar_load(c("debranched_seus", "clustree_tables", "large_clone_comparisons", "cluster_dictionary"))


# test1 <- find_all_diffex_from_clustree(clustree_tables[["SRR13884242"]], "SRR13884242", filtered_seus, clone_comparisons = large_clone_comparisons)

# debug(find_all_diffex_from_clustree)
debug(find_diffex_from_clustree)
# debug(find_diffex_bw_divergent_clusters)
# debug(make_cluster_comparison)

test1 <- find_all_diffex_from_clustree(clustree_tables[3], debranched_seus, clone_comparisons = large_clone_comparisons)

test0 <- imap(clustree_tables, find_all_diffex_from_clustree, filtered_seus, clone_comparisons = large_clone_comparisons)


# figures and tables ------------------------------

tar_load(c("table_all_diffex_clones", "oncoprint_enrich_clones_plots_gobp", "oncoprint_plots"))

# debug(copy_paper_data)

copy_paper_data(table_all_diffex_clones, oncoprint_enrich_clones_plots_gobp, oncoprint_plots)

# plot_seu_clusters_and_markers ------------------------------

tar_load(c("debranched_seus", "cluster_orders"))

# debug(plot_seu_clusters_and_markers)

plot_seu_clusters_and_markers(debranched_seus[1], cluster_orders)

# make_volcano_diffex_clones ------------------------------

tar_load(c("cis_diffex_clones_for_each_cluster", "cis_diffex_clones", "cis_diffex_clones_for_each_phase"))

debug(make_volcano_diffex_clones)

make_volcano_diffex_clones(
	cis_diffex_clones_for_each_cluster,
	"results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
	cis_diffex_clones,
	"results/diffex_bw_clones_large_in_segment.pdf"
)



# make_volcano_plots ------------------------------

res <- read_csv("results/numbat_sridhar/SRR13884243_cluster_clone_comparison_diffex_in_segment.csv")

myplot <- 
res %>% 
	dplyr::filter(clone_comparison == "2_v_1_1q+_16q-") %>% 
	dplyr::mutate(diffex_comparison = str_replace(str_extract(clone_comparison, "[0-9]_v_[0-9]"), "_v_", "_")) %>%
	dplyr::distinct(symbol, .keep_all = TRUE) %>% 
	tibble::column_to_rownames("symbol") %>%
	# dplyr::filter(location == "all") %>% 
	make_volcano_plots(sample_id = "asdf", mysubtitle = "asdf") %>%
	identity()
 
ggplot_build(myplot)$layout$panel_params[[1]]$y.range



# find_candidate_cis_in_clustree_diffexes ------------------------------

tar_load("clustree_diffexes")

debug(find_candidate_cis_in_clustree_diffexes)

find_candidate_cis_in_clustree_diffexes(clustree_diffexes)



# pull_cluster_orders------------------------------

tar_load("cluster_file")

debug(pull_cluster_orders)

pull_cluster_orders(cluster_file)

# pull_clustree_tables ------------------------------

tar_load("clustrees")

debug(pull_clustree_tables)

pull_clustree_tables(clustrees)




# enrichment_analysis ------------------------------

tar_load("oncoprint_plots")

diffex <- myreadxl(oncoprint_plots$cis$table)

# debug(enrichment_analysis)

enrichment_output <-
diffex$`2p+` %>%
  group_by(symbol) %>%
  dplyr::slice_max(abs_log2FC) %>%
  dplyr::distinct(symbol, .keep_all = TRUE) %>%
  tibble::column_to_rownames("symbol") %>%
  dplyr::select(-any_of(colnames(annotables::grch38))) %>%
  enrichment_analysis() %>%
  setReadable(
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  )

myplot <-
  enrichment_output %>%
  plot_enrichment()


# compare_markers ------------------------------

tar_load(c("filtered_seus", "cluster_orders", "subtype_markers"))

# undebug(compare_markers)

# compare_markers(filtered_seus[[8]], cluster_orders)

safe_compare_markers <- safely(compare_markers)

test1 <- map(filtered_seus, safe_compare_markers, cluster_orders) %>%
  # set_names(str_extract(filtered_seus, "SRR[0-9]*")) %>%
  # map("result") %>%
  # compact() %>%
  identity()

# names(test1) <- str_extract(filtered_seus, "SRR[0-9]*")

recurrence_table <-
  test1 %>%
  # map(dplyr::select, Gene.Name, clusters) %>%
  dplyr::bind_rows(.id = "sample_id") %>%
  dplyr::arrange(Gene.Name) %>%
  dplyr::group_by(Gene.Name) %>%
  dplyr::mutate(recurrence = n()) %>%
  dplyr::group_by(sample_id, Gene.Name) %>%
  dplyr::slice_max(abs(Average.Log.Fold.Change)) %>%
  dplyr::arrange(desc(recurrence), Gene.Name) %>%
  identity()

test3 <-
  recurrence_table %>%
  dplyr::select(sample_id, Gene.Name, Cluster) %>%
  tidyr::pivot_wider(names_from = "sample_id", values_from = "Cluster") %>%
  identity()

recurrence_table %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sample_id, clusters) %>%
  dplyr::mutate(Gene.Name = factor(Gene.Name, levels = unique(Gene.Name))) %>%
  ggplot(aes(x = sample_id, y = Gene.Name, fill = clusters)) +
  geom_tile()

ggsave("results/g1_cluster_overlap.pdf", width = 15, height = 50)


# debug(compare_cluster_continuous_var)

continuous_var = "subtype2"

# compare_cluster_continuous_var(filtered_seus[[3]], cluster_orders, continuous_var = continuous_var, gene_lists = subtype_markers)

test2 <- map(filtered_seus, compare_cluster_continuous_var, cluster_orders, continuous_var = continuous_var, gene_lists = subtype_markers)

names(test2) <- str_extract(filtered_seus, "SRR[0-9]*")

test5 <-
  test2  %>%
  dplyr::bind_rows(.id = "sample_id") %>%
  dplyr::group_by(sample_id, clusters) %>%
  dplyr::summarize(mean_cont = mean(!!sym(continuous_var))) %>%
  dplyr::arrange(desc(mean_cont)) %>%
  # dplyr::filter(str_detect(clusters, "g1_[0-9]")) %>%
  identity()

test5 %>%
  ggplot(aes(x = sample_id, y = mean_cont, color = clusters, label = clusters)) +
  geom_text()





# make_volcano_diffex_clones ------------------------------

tar_load(c(
  "cis_diffex_clones_for_each_cluster",
  "cis_diffex_clones",
  "cis_diffex_clones_for_each_phase"
))

# debug(make_volcano_diffex_clones)

make_volcano_diffex_clones(
  cis_diffex_clones_for_each_cluster,
  "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
  cis_diffex_clones,
  "results/diffex_bw_clones_large_in_segment.pdf",
  cis_diffex_clones_for_each_phase,
  "results/diffex_bw_clones_per_phase_large_in_segment.pdf"
)


# score_samples_for_celltype_enrichment ------------------------------

tar_load("celltype_markers")

# debug(score_samples_for_celltype_enrichment)
# debug(score_binary_celltype_markers)

score_samples_for_celltype_enrichment("output/seurat/SRR14800543_unfiltered_seu.rds", "output/seurat/SRR14800543_filtered_seu.rds", celltype_markers)

# heatmap_marker_genes ------------------------------

tar_load(c("filtered_common_genes", "subtype_markers"))

debug(heatmap_marker_genes)

heatmap_marker_genes("output/seurat/SRR13884242_filtered_seu.rds", filtered_common_genes, subtype_markers, "filtered_", marker_col = "SCT_snn_res.0.2", group.by = c("SCT_snn_res.0.2", "scna", "Phase"), col_arrangement = c("SCT_snn_res.0.2", "scna", "Phase"))



# pull_common_markers ------------------------------

tar_load(c("filtered_seus", "mps"))

debug(pull_common_markers)

pull_common_markers(filtered_seus, mps[["Cancer"]])

# score_and_heatmap_seu ------------------------------

tar_load("mps")

debug(score_and_heatmap_seu)

score_and_heatmap_seu("output/seurat/SRR13884242_filtered_seu.rds", mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv")

# score_and_vlnplot_seu ------------------------------

tar_load(c("subtype_markers"))

debug(score_and_vlnplot_seu)

score_and_vlnplot_seu("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds", subtype_markers, leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv", y_lim = 1.8, step = 0.2)


# pull_subtype_genes ------------------------------

debug(pull_subtype_genes)

pull_subtype_genes(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx")

# make_volcano_diffex_clones ------------------------------

# debug(make_volcano_diffex_clones)
# debug(annotate_cluster_membership)

tar_load(c("trans_diffex_clones_for_each_cluster", "trans_diffex_clones", "trans_diffex_clones_for_each_phase"))

make_volcano_diffex_clones(
  trans_diffex_clones_for_each_cluster,
  "results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf",
  trans_diffex_clones,
  "results/diffex_bw_clones_large_out_of_segment.pdf",
  trans_diffex_clones_for_each_phase,
  "results/diffex_bw_clones_per_phase_large_out_of_segment.pdf")

# inspect_oncoprints ------------------------------

tar_load("oncoprint_input_by_scna")

debug(inspect_oncoprints)

inspect_oncoprints(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans)


# plot_putative_marker_across_samples ------------------------------

tar_load(c("interesting_genes", "regressed_seus", "cluster_dictionary"))

# debug(plot_putative_marker_across_samples)

plot_putative_marker_across_samples(interesting_genes, unlist(regressed_seus), plot_type = FeaturePlot, group_by = "SCT_snn_res.0.4", cluster_dictionary, extension = "regressed")

# plot_effect_of_regression ------------------------------

# debug(plot_effect_of_regression)

plot_effect_of_regression("output/seurat/SRR13884242_filtered_seu.rds", "output/seurat/SRR13884242_regressed_seu.rds")

# tabulate_diffex_clones out_of_segment ------------------------------

debug(tabulate_diffex_clones)
# debug(annotate_cluster_membership)

tar_load(c("large_trans_diffex_clones_for_each_cluster", "large_trans_diffex_clones", "large_trans_diffex_clones_for_each_phase"))
tar_load(c("large_cis_diffex_clones_for_each_cluster", "large_cis_diffex_clones", "large_cis_diffex_clones_for_each_phase"))

tabulate_diffex_clones(large_cis_diffex_clones_for_each_cluster,
                       "results/diffex_bw_clones_per_cluster_large_in_segment.xlsx",
                       "results/diffex_bw_clones_per_cluster_large_in_segment_by_chr.xlsx",
                       large_cis_diffex_clones,
                       "results/diffex_bw_clones_large_in_segment.xlsx",
                       "results/diffex_bw_clones_large_in_segment_by_chr.xlsx",
                       large_cis_diffex_clones_for_each_phase,
                       "results/diffex_bw_clones_per_phase_large_in_segment.xlsx",
                       "results/diffex_bw_clones_per_phase_large_in_segment_by_chr.xlsx")

# score_and_vlnplot_seu ------------------------------

tar_load(c("subtype_markers"))

debug(score_and_vlnplot_seu)

score_and_vlnplot_seu("output/seurat/SRR13884242_filtered_seu.rds", subtype_markers, leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv")



# compile_subtype_violins ------------------------------

tar_load(c("interesting_samples", "subtype_violins"))

debug(compile_subtype_violins)

compile_subtype_violins(interesting_samples, subtype_violins)


# plot_celltype_predictions ------------------------------

tar_load("plae_ref")

debug(plot_celltype_predictions)

plot_celltype_predictions("output/seurat/SRR13884242_regressed_seu.rds", "SRR13884242", plae_ref)

# plot_plae_celltype_expression ------------------------------

clone_2_v_3_genes %>%
  plot_plae_celltype_expression()

c("STMN4", "GNG3", "MEG3", "CRABP2", "LY6H", "NDUFA4L2", "PDE6H",
  "PCM1", "ELAVL3", "ELAVL2") %>%
  plot_plae_celltype_expression()

plot_plae_celltype_expression(c("BCL11A", "ISOC1", "NRXN1", "CRYBG3", "RBP7", "PRDM1"))

plot_plae_celltype_expression(c("RORB", "CADPS", "FXYD3", "ARID1B"))

## clone distribution ------------------------------

filtered_seu_path <- "output/seurat/SRR13884242_filtered_seu.rds"
regressed_seu_path <- "output/seurat/SRR13884242_regressed_seu.rds"

regressed_seu <- readRDS(regressed_seu_path)

debug(table_distribution_of_clones_across_clusters)

table_distribution_of_clones_across_clusters(regressed_seu, "SRR13884242")

# groupGO ------------------------------

seu <- readRDS("output/seurat/SRR13884242_filtered_seu.rds")

seu <- find_all_markers(seu, "abbreviation")

test0 <-
seu@misc$markers$abbreviation$presto %>%
  dplyr::filter(`Adjusted.pvalue` < 0.05) %>%
  dplyr::arrange(`Average.Log.Fold.Change`) %>%
  dplyr::filter(`Average.Log.Fold.Change` > 0) %>%
  dplyr::left_join(annotables::grch38, by = c("Gene.Name" = "symbol")) %>%
  dplyr::mutate(entrez = as.character(entrez)) %>%
  dplyr::distinct(Cluster, entrez) %>%
  split(.$Cluster) %>%
  map(pull, entrez) %>%
  map(clusterProfiler::enrichGO, 'org.Hs.eg.db', ont="BP") %>%
  identity()

test1 <-
  test0 %>%
  map(as.data.frame) %>%
  map(arrange, desc(Count)) %>%
  map(head) %>%
  map(pull, Description) %>%
  identity()


# plot_effects_of_cell_cycle_regression -------------------------------

# plot_effects_of_cell_cycle_regression("")

# compile_cis_trans_enrichment_recurrence------------------------------

tar_load("oncoprint_enrich_clusters_hallmark")

debug(compile_cis_trans_enrichment_recurrence)

compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clusters_hallmark,
                                        cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_hallmark",
                                        trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_hallmark",
                                        cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_hallmark.xlsx",
                                        trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_hallmark.xlsx",
                                        by_cluster = TRUE)

# bulk_subtype_scoring ------------------------------

bulk_subtype_scoring <- function(seu_path, subtype_markers){
  sample_id = str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)

  pseudobulk <- GetAssayData(seu, assay = "gene", slot = "data") %>%
    rowSums()

  s1_score = summary(pseudobulk[names(pseudobulk) %in% subtype_markers$subtype1])

  s2_score = summary(pseudobulk[names(pseudobulk) %in% subtype_markers$subtype2])

  return(list("sample_id" = sample_id, "subtype1" = s1_score, "subtype2" = s2_score))

}

seu_paths <- dir_ls("output/seurat/", glob = "*filtered_seu.rds", recurse = TRUE)

test0 <- map_dfr(seu_paths, bulk_subtype_scoring, subtype_markers)

# compile_cis_trans_enrichment_recurrence ------------------------------

tar_load("oncoprint_enrich_clones_gobp")
#
# # debug("plot_enrichment_recurrence")
#
debug(compile_cis_trans_enrichment_recurrence)
debug(plot_enrichment_recurrence)

compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_gobp,
                                        cis_plot_file = "results/cis_enrichment_plots_by_clone_gobp",
                                        trans_plot_file = "results/trans_enrichment_plots_by_clone_gobp",
                                        cis_table_file = "results/cis_enrichment_tables_by_clone_gobp.xlsx",
                                        trans_table_file = "results/trans_enrichment_tables_by_clone_gobp.xlsx")


# plot violins ------------------------------

tar_load("subtype_markers")
debug(score_and_plot_seu)

score_and_plot_seu("output/seurat/SRR14800540_filtered_seu.rds", subtype_markers, leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv")



# 2 ------------------------------

tar_load(c("s1_markers", "s2_markers"))

test0 <-
  FetchData(seu, vars = s1_markers) %>%
  rowSums()


seu_paths <- dir_ls("output/seurat", glob = "*_filtered_seu.rds") %>%
  set_names(str_extract(., "SRR[0-9]*"))


seu <- Seurat::MetaFeature(seu, features = s2_markers, meta.name = "s2_score")

VlnPlot(seu, features = "s2_score", group.by = "abbreviation") +
  scale_y_log10()

seu <- Seurat::MetaFeature(seu, features = s1_markers, meta.name = "s1_score")

VlnPlot(seu, features = "s1_score", group.by = "clone_opt")

seu <- Seurat::MetaFeature(seu, features = c("TOP2A", "CENPF"), meta.name = "test0")

FeaturePlot(seu, features = "test0")

FeaturePlot(seu, features = "ESRRG")

