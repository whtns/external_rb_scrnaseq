# tar_target(
# plae_ref,
# generate_plae_ref()
# ),
# tar_target(
#   current_params,
#   retrieve_snakemake_params(numbat_rds_files),
#   pattern = map(numbat_rds_files),
#   iteration = "list"
# ),
#
# tar_target(
#   current_init_k,
#   retrieve_current_param(current_params, "init_k")
# ),
# tar_target(large_plot_files,
#   make_numbat_plot_files(numbat_rds_files, cluster_dictionary, clone_simplifications = large_clone_simplifications),
#   pattern = map(numbat_rds_files),
#   iteration = "list"
# ),
# tarchetypes::tar_file_fast(
# 	cluster_file,
# 	"data/raw_cluster_ids.csv"
# ),
  #
  # tar_target(merged_subtype_violins,
  #            score_and_vlnplot_seu("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds", subtype_markers, leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv", y_lim = 1.8, step = 0.2),
  # ),
  # tar_target(
  #   stemness_violins,
  #   score_and_vlnplot_seu(final_seus, stemness_markers[["smith"]]),
  #   pattern = map(final_seus),
  #   iteration = "list"
  # ),
    # tar_target(
  #   mp_violins,
  #   score_and_vlnplot_seu(final_seus, mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv"),
  #   pattern = map(final_seus),
  #   iteration = "list"
  # ),
    # tar_target(cc_plots_wo_arms,
  # 					 filter_cluster_save_seu(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, cells_to_remove, extension = "_filtered", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
  # 					 pattern = map(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions),
  # 					 iteration = "list"
  # ),
    # tar_target(debranched_seus,
  #   debranch_seus(final_seus, branch_dictionary)
  # ),
    # tar_target(filtered_scanpys,
  #            convert_seu_to_scanpy(final_seus),
  #            pattern = map(final_seus),
  #            iteration = "list"
  # ),
  # tar_target(regressed_scanpys,
  #            convert_seu_to_scanpy(regressed_seus),
  #            pattern = map(regressed_seus),
  #            iteration = "list"
  # ),
  # tar_target(clustify_plots,
  #            plot_celltype_predictions(regressed_seus, interesting_samples, plae_ref, group.by = "SCT_snn_res.0.6"),
  #            pattern = map(regressed_seus, interesting_samples),
  #            iteration = "list"
  # ),
    # tar_target(large_heatmaps,
  #   make_numbat_heatmaps_old(numbat_rds_files, large_filter_expressions, cluster_dictionary, p_min = 0.5, line_width = 0.1, extension = "_filtered"),
  #   pattern = map(numbat_rds_files, large_filter_expressions, cluster_dictionary),
  #   iteration = "list"
  # ),
    # tar_target(clustrees,
  # 					 make_clustrees_for_sample("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", mylabel = "sdfg", assay = "SCT")
  # ),
    # tar_target(cis_louvain_comparisons,
  #            find_diffex_bw_clones_for_each_louvain(final_seus, numbat_rds_files, large_clone_comparisons, location = "in_segment"),
  #            pattern = map(final_seus, numbat_rds_files, large_clone_comparisons),
  #            iteration = "list"
  # ),
  #
  # tar_target(trans_louvain_comparisons,
  #            find_diffex_bw_clones_for_each_louvain(final_seus, numbat_rds_files, large_clone_comparisons, location = "out_of_segment"),
  #            pattern = map(final_seus, numbat_rds_files, large_clone_comparisons),
  #            iteration = "list"
  # ),
  #
  # tar_target(all_louvain_comparisons,
  #            find_diffex_bw_clones_for_each_louvain(final_seus, numbat_rds_files, large_clone_comparisons, location = "all"),
  #            pattern = map(final_seus, numbat_rds_files, large_clone_comparisons),
  #            iteration = "list"
  # ),
    # # 16q
  # tar_target(integrated_seu_16q, seu_integrate_rbs(numbat_dir = "output/numbat_sridhar",
  #                                                  kept_samples = c('SRR14800534', 'SRR14800535', 'SRR14800536'),
  #                                                  cluster_dictionary)
  #            ),
  # tar_target(clone_comparisons_for_16q, list("2_v_1_16q-" = c("16b"), "3_v_2_1q+" = c("1b"))), #1q only GT 3; maybe reconsider
  # # 16q/1q
  # tar_target(integrated_seu_1q_16q, seu_integrate_rbs(numbat_dir = "output/numbat_sridhar",
  #                                                     kept_samples = c('SRR13884242', 'SRR13884243', 'SRR13884249', 'SRR14800534', 'SRR14800535', 'SRR14800536', 'SRR14800543'),
  #                                                     cluster_dictionary)
  #            ),
  # # 1q 6p
  # tar_target(integrated_seu_1q_6p, seu_integrate_rbs(numbat_dir = "output/numbat_sridhar",
  #                                                    kept_samples = c('SRR14800543', 'SRR17960481', 'SRR17960484'),
  #                                                    cluster_dictionary)
  #            ),
  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(cis_diffex_clones_for_each_cluster),
  #   pattern = map(cis_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(trans_diffex_clones_for_each_cluster),
  #   pattern = map(trans_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(cis_diffex_clones_for_each_cluster),
  #   pattern = map(cis_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(trans_diffex_clones_for_each_cluster),
  #   pattern = map(trans_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  # tar_target(large_diffex_bw_clusters_for_each_clone,
  #   find_diffex_bw_clusters_for_each_clone(final_seus, cluster_dictionary),
  #   pattern = map(final_seus, cluster_dictionary),
  #   iteration = "list"
  # ),
    # # check effect of regression
  # tar_target(heatmap_collages,
  # 					 plot_seu_marker_heatmap_all_resolutions(regressed_seus, numbat_rds_files, large_clone_simplifications),
  # 					 pattern = map(regressed_seus),
  # 					 iteration = "list"
  # ),
  # # cis diffex b/w clones whole tumor integrated_seus_1q
# tar_target(fig_07d,
# 					 plot_diffex_genes_on_split_integrated_seu()
# 					 ),
  # tar_target(large_numbat_bulk_clones_final, retrieve_numbat_plot_type(large_numbat_pdfs, "bulk_clones_final.pdf")),
  # tar_target(large_numbat_bulk_clones_final_pdf, qpdf::pdf_combine(large_numbat_bulk_clones_final, "results/numbat_sridhar_large/large_bulk_clones_final.pdf")),
  # # large_numbat_sample_pdfs ------------------------------
  # tar_target(large_numbat_sample_pdfs,
  #   reroute_done_to_results_pdf(numbat_rds_files, "_large"),
  #   pattern = map(numbat_rds_files),
  #   iteration = "list",
  #   format = "file"
  # ),
    # tar_target(clustree_diffex_test,
  # 					 find_all_diffex_from_clustree(clustree_tables[[1]], final_seus, clone_comparisons = large_clone_comparisons),
  # ),
  # tar_target(collage_compilation,
  #          qpdf::pdf_combine(heatmap_collages, "results/heatmap_collages.pdf")),
  # tar_target(large_clone_dist, retrieve_numbat_plot_type(large_plot_files, "exp_roll_clust.pdf")),
  # tar_target(large_numbat_expression_pdf, qpdf::pdf_combine(large_numbat_expression, "results/numbat_sridhar_large/large_numbat_expression.pdf")),
    # tar_target(
  #   oncoprint_input_by_scna_for_each_cluster_unfiltered,
  #   make_oncoprint_diffex_unfiltered(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, large_clone_comparisons, n_slice = 20)
  # ),
    # tar_target(
  #   num_diffex_cluster_scna_tally,
  #   tally_num_diffex(oncoprint_input_by_scna_for_each_cluster_unfiltered)
  # ),
    # clusters enrichment ------------------------------
  # # enrichment by cluster msigdb gobp
  # tar_target(
  #   oncoprint_enrich_clusters_gobp,
  #   enrich_oncoprints_clusters(large_filter_expressions,
  #     cluster_dictionary,
  #     debranched_ids,
  #     cis_diffex_clones_for_each_cluster,
  #     trans_diffex_clones_for_each_cluster,
  #     all_diffex_clones_for_each_cluster,
  #     large_clone_comparisons,
  #     gene_set = "gobp"
  #   )
  # ),
    # # enrichment plots and tables by cluster gobp
  # tar_target(
  #   oncoprint_enrich_clusters_plots_gobp,
  #   compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clusters_gobp,
  #     cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_gobp",
  #     trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_gobp",
  #     cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_gobp.xlsx",
  #     trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_gobp.xlsx",
  #     by_cluster = TRUE
  #   )
  # ),