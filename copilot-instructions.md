R is located at /opt/R/4.3.1/bin/R
Rscript is located at /opt/R/4.3.1/bin/Rscript
the targets store is _targets_r431
use the targets-pipeline skill
run the targets pipeline with /home/skevin/single_cell_projects/resources/external_rb_scrnaseq_proj/run_targets_bg.R

## Monitoring Targets Pipeline

When the user asks about targets pipeline progress, status, or to monitor the pipeline, use the monitoring script:
```bash
# Basic usage (uses default _targets_r431 store)
./scripts/monitor_targets.sh

# Specify different store
./scripts/monitor_targets.sh _targets

# Specify store and refresh interval (seconds)
./scripts/monitor_targets.sh _targets_r431 5
```

The script provides:
- Real-time progress percentage
- Count of completed, dispatched, and skipped targets
- List of currently active targets
- Error detection
- Recently completed targets with runtimes

For one-time status checks without continuous monitoring, use:
```bash
Rscript -e 'library(targets); tar_config_set(store = "_targets_r431"); progress <- tar_progress(); table(progress$progress)'
```
## Function Index

| Function | File |
|---|---|
| add_batch_hash_metadata | R/metadata_functions_1.R |
| add_cluster_label | R/plot_functions_26.R |
| add_hash_metadata | R/metadata_functions_1.R |
| add_hypoxia_score | R/heatmap_functions.R |
| add_image_slide | R/plot_functions_10.R |
| add_plots_slide | R/plot_functions_10.R |
| all_same_sign | R/plot_functions_8.R |
| annotables_to_grange | R/numbat_functions_10.R |
| annotate_cell_cycle_singleR | R/plot_functions_20.R |
| annotate_cell_cycle_without_1q | R/numbat_functions_12.R |
| annotate_cluster_membership | R/diffex_functions_10.R |
| annotate_cluster_scna_percentage | R/plot_functions_12.R |
| annotate_percent_segment_diffex | R/plot_functions_32.R |
| annotate_seu_with_rb_subtype_gene_expression | R/plot_functions_21.R |
| append_clone_nums | R/plot_functions_38.R |
| arrange_by_recurrence | R/plot_functions_49.R |
| assign_designated_phase_clusters | R/plot_functions_8.R |
| assign_phase_cluster_at_resolution | R/plot_functions_1.R |
| assign_phase_clusters | R/plot_functions_50.R |
| binom_test_clone_percent | R/numbat_functions_14.R |
| browse_celltype_expression | R/plot_functions_29.R |
| bulk_extract_seu_metadata | R/seu_metadata_db.R |
| calc_silhouette | R/plot_functions_52.R |
| calculate_6p_percent_affected | R/numbat_functions_14.R |
| calculate_arm_percent | R/numbat_functions_14.R |
| calculate_avg_cluster_distance | R/plot_functions_10.R |
| calculate_clone_distribution | R/plot_functions_1.R |
| check_cluster_marker_gene | R/numbat_functions_13.R |
| check_integrated_cluster_numbers | R/heatmap_functions.R |
| check_merged_metadata | R/plot_functions_8.R |
| check_MYCN_mean_expression | R/plot_functions_20.R |
| check_table_set | R/plot_functions_48.R |
| chi_sq_daughter_clusters | R/plot_functions_46.R |
| clean_diffex | R/plot_functions_51.R |
| clone_cc_plots_by_scna | R/plot_functions_10.R |
| clone_diff_per_cluster | R/plot_functions_25.R |
| cluster_columns | R/numbat_functions_12.R |
| cluster_diff_per_clone | R/plot_functions_25.R |
| collate_clone_distribution_tables | R/plot_functions_47.R |
| collect_all_markers | R/plot_functions_29.R |
| collect_clusters_from_seus | R/plot_functions_29.R |
| collect_markers | R/plot_functions_29.R |
| collect_study_metadata | R/diffex_functions_10.R |
| color_clustree_by_clone | R/plot_functions_45.R |
| compare_cluster_composition | R/seu_metadata_db.R |
| compare_cluster_continuous_var | R/plot_functions_29.R |
| compare_corresponding_enrichments | R/plot_functions_11.R |
| compare_diffex_cluster_by_clone | R/plot_functions_25.R |
| compare_enrichment | R/plot_functions_31.R |
| compare_infercnv | R/enrichment_functions_2.R |
| compare_markers | R/plot_functions_29.R |
| compare_per_cluster_and_total_clone_diffex | R/diffex_functions_2.R |
| compile_cis_trans_enrichment_recurrence | R/plot_functions_36.R |
| compile_cis_trans_enrichment_recurrence_by_cluster | R/plot_functions_36.R |
| compile_subtype_violins | R/plot_functions_29.R |
| compplot_feature_and_clusters | R/plot_functions_7.R |
| convert_numbat_pngs | R/enrichment_functions_2.R |
| convert_seu_to_scanpy | R/plot_functions_42.R |
| convert_v5_to_v3 | R/numbat_functions_12.R |
| convert_volcano_to_plotly | R/plot_functions_33.R |
| debranch_seus | R/numbat_functions_1.R |
| derive_pseudobulk_subtype_scores | R/plot_functions_39.R |
| diffex_by_cluster | R/enrichment_functions_1.R |
| diffex_cells | R/diffex_functions_1.R |
| diffex_groups | R/diffex_functions_1.R |
| diffex_groups_old | R/plot_functions_23.R |
| diffex_per_cluster | R/plot_functions_18.R |
| dotplot_diffex | R/plot_functions_14.R |
| dotplot_recurrent_genes | R/plot_functions_14.R |
| drop_bad_cells | R/plot_functions_36.R |
| drop_cc_genes | R/enrichment_functions_2.R |
| drop_cluster_orders_table | R/metadata_functions_1.R |
| drop_empty_plots | R/plot_functions_26.R |
| drop_mt_cluster | R/plot_functions_14.R |
| drop_y_axis | R/plot_functions_12.R |
| dummy_cluster_order | R/plot_functions_19.R |
| encode_cluster_order_to_hash_table | R/metadata_functions_1.R |
| enrich_by_cluster | R/enrichment_functions_2.R |
| enrich_diffex_by_cluster | R/enrichment_functions_2.R |
| enrichment_analysis | R/enrichment_functions_1.R |
| enrich_oncoprints | R/plot_functions_35.R |
| enrich_oncoprints_clusters | R/plot_functions_35.R |
| extract_full_segmentation | R/plot_functions_19.R |
| extract_seu_metadata | R/seu_metadata_db.R |
| filter_cluster_save_seu | R/plot_functions_3.R |
| filter_diffex | R/diffex_functions_3.R |
| filter_enrichment_result | R/plot_functions_11.R |
| filter_input_by_region | R/diffex_functions_3.R |
| filter_input_by_scna | R/plot_functions_33.R |
| filter_numbat_cells | R/plot_functions_7.R |
| filter_oncoprint_diffex | R/diffex_functions_3.R |
| filter_ora_result_by_terms | R/enrichment_functions_1.R |
| filter_phylo_plot | R/plot_functions_5.R |
| filter_sample_qc | R/plot_functions_3.R |
| filter_seu_to_cluster | R/enrichment_functions_1.R |
| find_all_diffex_from_clustree | R/plot_functions_48.R |
| find_candidate_all_in_clustree_diffexes | R/plot_functions_49.R |
| find_candidate_cis_in_clustree_diffexes | R/plot_functions_49.R |
| find_candidate_trans_in_clustree_diffexes | R/plot_functions_49.R |
| find_cc_genes_by_arm | R/plot_functions_52.R |
| find_cluster_pairwise_distance | R/plot_functions_11.R |
| find_diffex_bw_clones_for_each_cluster | R/plot_functions_26.R |
| find_diffex_bw_clones_for_each_cluster_integrated | R/plot_functions_36.R |
| find_diffex_bw_clusters_for_each_clone | R/plot_functions_25.R |
| find_diffex_bw_divergent_clusters | R/plot_functions_36.R |
| find_diffex_clones | R/plot_functions_30.R |
| find_diffex_clones_between_corresponding_states | R/plot_functions_31.R |
| find_diffex_clones_in_phase | R/plot_functions_30.R |
| find_diffex_clones_integrated | R/plot_functions_31.R |
| find_diffex_clusters_between_corresponding_states | R/plot_functions_31.R |
| find_diffex_from_clustree | R/plot_functions_48.R |
| find_genes_by_arm | R/plot_functions_52.R |
| generate_plae_ref | R/plot_functions_43.R |
| get_arms_ranges | R/numbat_functions_10.R |
| get_cell_qc_values | R/seu_metadata_db.R |
| get_chr_arms | R/numbat_functions_10.R |
| get_cluster_composition | R/seu_metadata_db.R |
| get_cluster_markers | R/seu_metadata_db.R |
| get_merged_metadata | R/plot_functions_23.R |
| get_metadata_summary | R/seu_metadata_db.R |
| get_qc_metrics | R/seu_metadata_db.R |
| gse_plot_from_clone_diffex | R/plot_functions_27.R |
| gse_plot_from_cluster_diffex | R/plot_functions_26.R |
| heatmap_from_cluster_diffex | R/plot_functions_17.R |
| heatmap_marker_genes | R/plot_functions_36.R |
| init_seu_metadata_db | R/seu_metadata_db.R |
| inspect_oncoprints | R/plot_functions_34.R |
| integration_by_scna_clones | R/integration_functions.R |
| interleave_lists | R/plot_functions_8.R |
| label_enrichment_by_cis | R/plot_functions_19.R |
| list_seu_objects | R/seu_metadata_db.R |
| load_and_save_hypoxia_score | R/heatmap_functions.R |
| make_all_numbat_plots | R/plot_functions_36.R |
| make_annoHighlight_from_consensus | R/plot_functions_16.R |
| make_cc_plot | R/plot_functions_9.R |
| make_cell_cycle_plot | R/enrichment_functions_1.R |
| make_clone_comparison | R/diffex_functions_2.R |
| make_clone_comparison_integrated | R/diffex_functions_2.R |
| make_clone_distribution_figure | R/plot_functions_26.R |
| make_cluster_comparison | R/diffex_functions_2.R |
| make_cluster_comparisons_by_phase_for_disctinct_clones | R/plot_functions_51.R |
| make_clustree_for_clone_comparison | R/plot_functions_15.R |
| make_clustrees_for_sample | R/plot_functions_20.R |
| make_corresponding_states_dictionary | R/plot_functions_16.R |
| make_dotplot | R/enrichment_functions_1.R |
| make_expression_heatmap_comparison | R/plot_functions_28.R |
| make_faded_umap_plots | R/plot_functions_46.R |
| make_filepaths_unique_in_hashes_table | R/metadata_functions_1.R |
| make_filtered_numbat_plots | R/plot_functions_6.R |
| make_gistic_plot | R/plot_functions_19.R |
| make_integrated_collage | R/plot_functions_20.R |
| make_integrated_numbat_plots | R/plot_functions_15.R |
| make_leng_sce | R/plot_functions_20.R |
| make_numbat_heatmaps | R/plot_functions_4.R |
| make_numbat_plot_files | R/plot_functions_22.R |
| make_oncoprint_diffex | R/diffex_functions_3.R |
| make_oncoprint_diffex_unfiltered | R/plot_functions_34.R |
| make_oncoprint_plots | R/plot_functions_34.R |
| make_pairwise_plots | R/plot_functions_44.R |
| make_pdf_montages | R/plot_functions_28.R |
| make_pearls_plot | R/plot_functions_21.R |
| make_rb_scna_ideograms | R/plot_functions_16.R |
| make_signed_dotplot | R/enrichment_functions_1.R |
| make_table_s01 | R/plot_functions_12.R |
| make_table_s02 | R/numbat_functions_13.R |
| make_table_s03 | R/numbat_functions_13.R |
| make_table_s05 | R/numbat_functions_10.R |
| make_table_s07 | R/plot_functions_19.R |
| make_table_s08 | R/numbat_functions_13.R |
| make_volcano_diffex_clones | R/plot_functions_32.R |
| make_volcano_plots | R/plot_functions_16.R |
| merge_orders | R/plot_functions_51.R |
| montage_images | R/plot_functions_52.R |
| myreadxl | R/plot_functions_50.R |
| numeric_col_fun | R/plot_functions_19.R |
| ora_analysis | R/enrichment_functions_1.R |
| ora_effect_of_regression | R/plot_functions_41.R |
| paired_plots_to_pptx | R/plot_functions_10.R |
| parse_metadata_summary | R/seu_metadata_db.R |
| pick_max_cluster | R/plot_functions_41.R |
| plot_cc | R/plot_functions_52.R |
| plot_cc_space_plot | R/plot_functions_53.R |
| plot_celltype_predictions | R/plot_functions_43.R |
| plot_clone_cc_plots | R/plot_functions_52.R |
| plot_clone_pearls | R/plot_functions_52.R |
| plot_clone_tree | R/plot_functions_48.R |
| plot_clone_tree_from_path | R/plot_functions_47.R |
| plot_cluster_markers_by_cell_type | R/plot_functions_7.R |
| plot_clustree_per_comparison | R/plot_functions_45.R |
| plot_corresponding_clusters_diffex_heatmaps | R/plot_functions_9.R |
| plot_corresponding_clusters_diffex_volcanos | R/plot_functions_9.R |
| plot_corresponding_enrichment | R/plot_functions_11.R |
| plot_diffex_genes_on_split_integrated_seu | R/plot_functions_15.R |
| plot_distribution_of_clones_across_clusters | R/plot_functions_2.R |
| plot_distribution_of_clones_pearls | R/plot_functions_21.R |
| plot_effect_of_filtering | R/plot_functions_40.R |
| plot_effect_of_regression | R/plot_functions_40.R |
| plot_enrichment | R/enrichment_functions_1.R |
| plot_enrichment_per_scna | R/plot_functions_38.R |
| plot_enrichment_recurrence | R/plot_functions_38.R |
| plot_enrichment_recurrence_by_cluster | R/plot_functions_38.R |
| plot_feature_across_seus | R/plot_functions_5.R |
| plot_feature_in_seu | R/plot_functions_52.R |
| plot_fig_04_afterall | R/figure_functions.R |
| plot_fig_03 | R/figure_functions.R |
| plot_fig_03_09 | R/plot_functions_11.R |
| plot_fig_03_afterall | R/figure_functions.R |
| plot_fig_04_05 | R/figure_functions.R |
| plot_fig_04_05_panels | R/figure_functions.R |
| plot_fig_04_07 | R/figure_functions.R |
| plot_fig_07_08 | R/figure_functions.R |
| plot_fig_08d | R/plot_functions_15.R |
| plot_fig_09_10 | R/figure_functions.R |
| plot_fig_s06a | R/figure_functions.R |
| plot_fig_s08 | R/figure_functions.R |
| plot_fig_s09 | R/figure_functions.R |
| plot_fig_s10 | R/figure_functions.R |
| plot_fig_s20 | R/figure_functions.R |
| plot_fig_s25 | R/plot_functions_20.R |
| plot_figure_collage | R/figure_functions.R |
| plot_filtering_timeline | R/plot_functions_44.R |
| plot_gene_clone_trend | R/plot_functions_30.R |
| plot_hypoxia_score | R/plot_functions_19.R |
| plot_integrated_1q_fig | R/plot_functions_12.R |
| plot_integrated_1q_fig_low_hypoxia | R/plot_functions_12.R |
| plot_markers_featureplot | R/plot_functions_36.R |
| plot_markers_in_sample | R/plot_functions_24.R |
| plot_merged_heatmap | R/plot_functions_9.R |
| plot_mt_v_nUMI | R/plot_functions_30.R |
| plot_numbat | R/plot_functions_1.R |
| plot_numbat_w_phylo | R/plot_functions_1.R |
| plot_pcnv_by_nsnp | R/diffex_functions_1.R |
| plot_pcnv_by_reads | R/enrichment_functions_2.R |
| plot_phase_distribution_by_scna | R/plot_functions_43.R |
| plot_phase_distribution_of_all_samples_by_scna | R/plot_functions_30.R |
| plot_phase_wo_arm | R/plot_functions_53.R |
| plot_plae_celltype_expression | R/plot_functions_30.R |
| plot_putative_marker_across_samples | R/plot_functions_25.R |
| plot_rb_subtype_expression | R/plot_functions_7.R |
| plot_recurrence | R/plot_functions_34.R |
| plot_scna_violins | R/plot_functions_36.R |
| plot_seu_clusters_and_markers | R/figure_functions.R |
| plot_seu_gene_heatmap | R/plot_functions_8.R |
| plot_seu_marker_heatmap | R/plot_functions_19.R |
| plot_seu_marker_heatmap_all_resolutions | R/plot_functions_46.R |
| plot_seu_marker_heatmap_by_scna | R/plot_functions_19.R |
| plot_seu_marker_heatmap_by_scna_ara | R/plot_functions_6.R |
| plot_seu_marker_heatmap_integrated | R/plot_functions_6.R |
| plot_several_diffex_clones_in_phase | R/plot_functions_6.R |
| plot_study_cell_stats | R/plot_functions_6.R |
| plot_study_metadata | R/plot_functions_39.R |
| plot_tcga_gistic | R/plot_functions_19.R |
| plot_variability_at_SCNA | R/plot_functions_1.R |
| plot_venn_w_genes | R/figure_functions.R |
| poster_plot_markers | R/plot_functions_52.R |
| prep_for_enrichment | R/plot_functions_38.R |
| prep_seu_branch | R/numbat_functions_1.R |
| prep_unfiltered_seu | R/plot_functions_3.R |
| print_table_tally | R/plot_functions_52.R |
| pub_violin | R/plot_functions_36.R |
| pull_assay_data | R/plot_functions_39.R |
| pull_branches | R/numbat_functions_1.R |
| pull_cells_matching_expression | R/plot_functions_29.R |
| pull_cluster_comparisons | R/plot_functions_51.R |
| pull_cluster_markers | R/plot_functions_44.R |
| pull_cluster_orders | R/plot_functions_46.R |
| pull_clustree_tables | R/numbat_functions_1.R |
| pull_common_markers | R/plot_functions_44.R |
| pull_hallmark_genes | R/plot_functions_21.R |
| pull_scna_segments | R/plot_functions_18.R |
| pull_segment_from_nb | R/diffex_functions_2.R |
| pull_segments | R/diffex_functions_2.R |
| pull_stem_cell_markers | R/plot_functions_37.R |
| pull_subtype_genes | R/plot_functions_37.R |
| rb_scna_frequency_in_tcga_by_cancer_type | R/plot_functions_17.R |
| read_all_hypoxia_scores | R/plot_functions_3.R |
| read_batch_hashes | R/metadata_functions_1.R |
| read_cells_to_remove | R/plot_functions_28.R |
| read_cluster_dictionary | R/plot_functions_1.R |
| read_cluster_orders_table | R/metadata_functions_1.R |
| read_giotti_genes | R/plot_functions_43.R |
| read_liu_lu_supp_tables | R/plot_functions_50.R |
| read_mps | R/diffex_functions_10.R |
| read_resolution_dictionary | R/metadata_functions_1.R |
| read_seu_hash | R/metadata_functions_1.R |
| read_seu_path | R/metadata_functions_1.R |
| read_zinovyev_genes | R/plot_functions_43.R |
| reference_plae_celltypes | R/plot_functions_42.R |
| regress_filtered_seu | R/plot_functions_3.R |
| remove_non_tumor_cells | R/plot_functions_33.R |
| reorder_within | R/plot_functions_17.R |
| reroute_done_to_results_pdf | R/numbat_utils_1.R |
| rescale_and_clean_plots | R/plot_functions_37.R |
| retrieve_cell_stats | R/diffex_functions_10.R |
| retrieve_current_param | R/metadata_functions_1.R |
| retrieve_genes_in_cis | R/numbat_functions_13.R |
| retrieve_numbat_plot_type | R/numbat_utils_1.R |
| retrieve_numbat_rds_files | R/numbat_utils_1.R |
| retrieve_numbat_seurat | R/plot_functions_6.R |
| retrieve_segmentation | R/numbat_functions_14.R |
| retrieve_seus | R/numbat_utils_1.R |
| retrieve_snakemake_params | R/metadata_functions_1.R |
| run_enrich_go | R/enrichment_functions_1.R |
| run_speckle | R/plot_functions_45.R |
| run_speckle_for_set | R/plot_functions_45.R |
| run_velocity | R/plot_functions_50.R |
| safe_plot_numbat | R/plot_functions_1.R |
| safe_plot_numbat_w_phylo | R/plot_functions_1.R |
| sample_specific_analyses_of_tumors_with_scna_subclones_after_integration | R/figure_functions.R |
| save_cc_space_plot_from_path | R/plot_functions_48.R |
| save_clone_tree_from_path | R/plot_functions_47.R |
| scale_x_reordered | R/plot_functions_17.R |
| scale_y_reordered | R/plot_functions_17.R |
| score_and_heatmap_seu | R/plot_functions_37.R |
| score_and_vlnplot_seu | R/plot_functions_36.R |
| score_binary_celltype_clusters | R/plot_functions_41.R |
| score_binary_celltype_markers | R/plot_functions_41.R |
| score_chrom_instability | R/plot_functions_42.R |
| score_clusters_up_down | R/plot_functions_47.R |
| score_filtration | R/plot_functions_29.R |
| score_pca | R/plot_functions_41.R |
| score_samples_for_celltype_enrichment | R/plot_functions_41.R |
| score_samples_for_rod_enrichment | R/plot_functions_39.R |
| score_stachelek | R/plot_functions_42.R |
| score_whole_pseudobulks | R/plot_functions_39.R |
| select_16q_clones | R/heatmap_functions.R |
| select_1q_clones | R/heatmap_functions.R |
| select_genes_from_arbitrary_diffex | R/plot_functions_49.R |
| select_genes_to_plot | R/plot_functions_12.R |
| select_scna_diffex | R/diffex_functions_3.R |
| send_pipeline_notification | R/send_notification.R |
| set_final_seus | R/plot_functions_50.R |
| seu_complex_heatmap | R/plot_functions_53.R |
| seu_complex_heatmap2 | R/plot_functions_53.R |
| seu_factor_heatmap | R/numbat_functions_12.R |
| seu_gene_heatmap | R/plot_functions_53.R |
| seu_integrate_rbs | R/plot_functions_41.R |
| shifter | R/plot_functions_21.R |
| simplify_gt | R/plot_functions_41.R |
| simplify_gt_col | R/plot_functions_15.R |
| split_by_cnv | R/diffex_functions_1.R |
| split_label_line | R/plot_functions_41.R |
| split_seu_by_branch | R/numbat_functions_1.R |
| split_vector_into_overlapping_chunks | R/plot_functions_20.R |
| subset_seu_by_clones | R/integration_functions.R |
| subset_seu_by_expression | R/plot_functions_19.R |
| subset_to_16q | R/heatmap_functions.R |
| subset_to_1q | R/heatmap_functions.R |
| subset_to_high_hypoxia | R/heatmap_functions.R |
| subset_to_low_hypoxia | R/heatmap_functions.R |
| summarize_clusters | R/plot_functions_47.R |
| table_and_plot_enrichment | R/plot_functions_51.R |
| table_cluster_markers | R/plot_functions_22.R |
| table_distribution_of_clones_across_clusters | R/plot_functions_22.R |
| tabulate_clone_comparisons | R/plot_functions_33.R |
| tabulate_clone_percent | R/numbat_functions_14.R |
| tabulate_clone_percent_integrated | R/numbat_functions_10.R |
| tabulate_diffex_clones | R/plot_functions_32.R |
| tally_kooi_candidates | R/diffex_functions_10.R |
| tally_num_diffex | R/diffex_functions_10.R |
| tidy_eval_arrange | R/plot_functions_45.R |
| upsert_resolution_dictionary | R/metadata_functions_1.R |
| volcano_plot_clone_clusters | R/plot_functions_32.R |
| volcano_plot_clones | R/plot_functions_32.R |
| window | R/plot_functions_26.R |
| wrap_scna_labels | R/plot_functions_44.R |

 why do i keep getting this error while running the r targets package? arguments imply differing number of rows: 1, 0
Gemini said

This error is a classic R headache, but in the context of the targets package, it usually signals a mismatch between what a function expects to receive and what it actually gets during the pipeline execution.

The error arguments imply differing number of rows: 1, 0 specifically means you are trying to combine or create a data frame (or tibble) where one column has 1 value and another is completely empty (0 rows).
Common Causes in targets
1. The "Empty Result" Trap

This is the most frequent culprit. One of your upstream targets is likely returning an empty vector or a data frame with zero rows (perhaps due to a filter that didn't match anything). When a downstream function tries to column-bind that empty result with a single value (like a label or a constant), R throws this error.
2. Issues with tar_map or Branching

If you are using dynamic branching or tar_map, the error often occurs when:

    The input data for the map is empty.

    The grid of parameters you provided has a null or empty value for one of the iterations.

3. Unexpected NULL values

If a function in your pipeline fails silently and returns NULL, and you later try to put that NULL into a data.frame(), R interprets the NULL as having 0 rows.
How to Debug It

Since targets runs in a clean environment, debugging can be tricky. Try these steps:

    Inspect the Upstream Target:
    Run tar_read(your_upstream_target_name) in your console to see exactly what the data looks like right before the error happens. Check if it's empty (nrow == 0).

    Use tar_meta():
    Check the metadata for the failing target to see if there are warnings or if the "bytes" size is suspiciously small.

    The "Workspace" Method:
    targets can save the exact state of a failure. Add tar_option_set(workspace_on_error = TRUE) to your _targets.R file. After it fails, run tar_workspace(name_of_failed_target) to load the exact environment and run the code line-by-line.

Example of the Bug
R

# This causes the error:
data.frame(id = "A", value = numeric(0)) 
# Error: arguments imply differing number of rows: 1, 0

Potential Fixes

    Add an if check: Inside your functions, check if the data is empty before trying to build a data frame.

    Return NA instead of empty: If a calculation finds nothing, returning a single NA (which has a length of 1) will often prevent the row mismatch.

    Validate Inputs: Use stopifnot() or the checkmate package at the start of your functions to catch empty inputs before they hit the data frame stage.