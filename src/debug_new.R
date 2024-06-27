
source("packages.R")
source("functions.R")
library(targets)


# plot_seu_marker_heatmap ------------------------------

tar_load(c("debranched_seus_6p", "debranched_seus_2p", "debranched_seus_16q", "debranched_seus_1q", "cluster_orders", "large_clone_simplifications", "numbat_rds_files", "rb_scna_samples", "large_clone_comparisons"))

# debug(plot_seu_marker_heatmap_by_scna)
# debug(plot_distribution_of_clones_across_clusters)

max_test0 <-  plot_seu_marker_heatmap_by_scna(debranched_seus_2p[[4]], cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", min_cells_per_cluster = 10)

min_test0 <- plot_seu_marker_heatmap_by_scna(debranched_seus_1q[[3]], cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 10, label = "minimum_resolution")

test0 <- plot_seu_marker_heatmap_by_scna(debranched_seus_1q[[3]], cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 10, label = "medium_resolution")



# plot_seu_marker_heatmap_all_resolutions ------------------------------

# debug(plot_seu_marker_heatmap_all_resolutions)
# debug(plot_seu_marker_heatmap)
# debug(plot_clone_tree)

tar_load(c("debranched_seus", "debranched_cluster_orders", "large_clone_simplifications", "numbat_rds_files"))

test0 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[10]], numbat_rds_files, large_clone_simplifications, debranched_cluster_orders)

test1 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[11]], numbat_rds_files, large_clone_simplifications, debranched_cluster_orders)

test2 <- plot_seu_marker_heatmap_all_resolutions(debranched_seus[[12]], numbat_rds_files, large_clone_simplifications, debranched_cluster_orders)



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


# make_clustree_for_speckle_set ------------------------------

tar_load(c("debranched_seus", "debranched_samples"))

debranched_seus <-
	debranched_seus %>%
	set_names(str_extract(debranched_seus, "SRR[0-9]*")) %>%
	identity()

debug(make_clustrees_for_sample)
debug(make_clustree_for_clone_comparison)
# debug(plot_clustree_per_comparison)
debug(color_clustree_by_clone)

make_clustrees_for_sample(debranched_seus[[3]], mylabel = debranched_samples[[3]], assay = "SCT")


# plot_clone_pearls ------------------------------

tar_load(c("final_seus", "cluster_orders"))

debug(plot_clone_pearls)
debug(plot_distribution_of_clones_pearls)

plot_clone_pearls(final_seus[3], cluster_orders[3])



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


# find_diffex_bw_clones_for_each_cluster ------------------------------

# debug(find_diffex_bw_clones_for_each_cluster)
# debug(clone_diff_per_cluster)

tar_load(c("debranched_seus", "numbat_rds_files",  "large_clone_comparisons", "debranched_cluster_orders"))

find_diffex_bw_clones_for_each_cluster(debranched_seus[[1]], numbat_rds_files,  large_clone_comparisons, debranched_cluster_orders, location = "out_of_segment")



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





# make oncoprint plots ------------------------------

tar_load(c("oncoprint_input_by_scna"))

debug(make_oncoprint_plots)

debug(plot_recurrence)

make_oncoprint_plots(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans, oncoprint_input_by_scna$all)


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

# plot_effect_of_filtering_regression ------------------------------

# tar_load("plot_effect_of_filtering_regression")

debug(plot_effect_of_filtering)

plot_effect_of_filtering("output/seurat/SRR14800541_unfiltered_seu.rds", "output/seurat/SRR14800541_filtered_seu.rds")





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

# ora_effect_of_regression ------------------------------

debug(ora_effect_of_regression)

ora_effect_of_regression("output/seurat/SRR13884242_filtered_seu.rds", "output/seurat/SRR13884242_regressed_seu.rds")

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

# effect of regression ------------------------------

debug(plot_effect_of_regression)

plot_effect_of_regression("output/seurat/SRR14800540_filtered_seu.rds",
                          "output/seurat/SRR14800540_regressed_seu.rds")


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

