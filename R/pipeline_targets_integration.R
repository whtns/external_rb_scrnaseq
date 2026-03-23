# Integration analysis targets: corresponding-state diffex (2p, 6p, combined),
# SCNA collages, oncoprint targets, enrichment, pseudobulks, and figures 09/10.
# Defines: pipeline_targets_integration (spliced into tar_plan in _targets.R)

pipeline_targets_integration <- list(

  # --- collages: sample-specific SCNA heatmap collages ---

  tar_target(collages_2p,
    plot_seu_marker_heatmap_by_scna(
      unlist(debranched_seus_2p), cluster_orders, numbat_rds_files, large_clone_simplifications,
      rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p"
    ),
    pattern = map(debranched_seus_2p),
    iteration = "list"
  ),

  # Sample-specific analyses of tumors with 2p+ subclones without integration (Fig. S4.9)
  tar_target(fig_s04_09,
    qpdf::pdf_combine(collages_2p, "results/fig_s04_09.pdf")
  ),

  tar_target(collages_6p,
    plot_seu_marker_heatmap_by_scna(
      unlist(debranched_seus_6p), cluster_orders, numbat_rds_files, large_clone_simplifications,
      rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p"
    ),
    pattern = map(debranched_seus_6p),
    iteration = "list"
  ),

  # Sample-specific analyses of tumors with 16q- subclones without integration.
  tar_target(collages_16q,
    plot_seu_marker_heatmap_by_scna(
      unlist(debranched_seus_16q), cluster_orders, numbat_rds_files, large_clone_simplifications,
      rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q"
    ),
    pattern = map(debranched_seus_16q),
    iteration = "list"
  ),

  tar_target(collected_scna_collages,
    qpdf::pdf_combine(
      unlist(list(collages_2p, collages_6p, collages_16q)),
      "results/fig_s17.pdf")
  ),

  # --- corresponding states: shared seu paths ---

  tar_target(corresponding_seus_2p,
    c(
      # "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds",
      # "SRR13884247_branch_6_filtered_seu" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
      "SRR13884248_filtered_seu_2p" = "output/seurat/SRR13884248_filtered_seu_2p.rds",
      "SRR17960484_filtered_seu_2p.rds" = "output/seurat/SRR17960484_filtered_seu_2p.rds"
      # "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"
    )
  ),

  tar_target(corresponding_seus_6p,
    c(
      "SRR13884248_filtered_seu_6p.rds" = "output/seurat/SRR13884248_filtered_seu_6p.rds",
      "SRR17960484_filtered_seu_6p.rds" = "output/seurat/SRR17960484_filtered_seu_6p.rds"
    )
  ),

  tar_target(corresponding_seus,
    list(
      # "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds",
      "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
      "output/seurat/SRR13884248_filtered_seu_2p.rds",
      "output/seurat/SRR17960484_filtered_seu_2p.rds",
      "output/seurat/SRR13884247_filtered_seu.rds",
      "output/seurat/SRR17960484_filtered_seu_6p.rds",
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds",
      "output/seurat/integrated_6p/integrated_seu_6p_duo.rds",
      "output/seurat/SRR13884248_filtered_seu_6p.rds"
    )
  ),

  tar_target(corresponding_states_dictionary,
    make_corresponding_states_dictionary(tibble::tribble(
      ~file_name, ~w_scna, ~wo_scna, ~scna_of_interest,
      # # "SRR13884246_branch_5_filtered_seu_2p.rds", "g1_1-g1_0-g1_3-g1_11", "g1_2-g1_5", "2p",
      # # "SRR13884246_branch_5_filtered_seu_2p.rds",     "s_g2_7-s_9",      "s_g2_8", "2p",

      "SRR13884247_branch_6_filtered_seu.rds", "g1_3-g1_6", "g1_5", "2p",
      # "SRR13884247_branch_6_filtered_seu.rds", "g1_3", "g1_5", "2p",
      # "SRR13884247_branch_6_filtered_seu.rds", "g1_6", "g1_5", "2p",

      "SRR13884248_filtered_seu_2p.rds",   "g1_1",    "g1_6-g1_7", "2p",
      # "SRR13884248_filtered_seu_2p.rds",   "g1_1",    "g1_6", "2p",
      # "SRR13884248_filtered_seu_2p.rds",   "g1_1",    "g1_7", "2p",

      "SRR17960484_filtered_seu_2p.rds",     "g1_3",    "g1_1-g1_5", "2p",
      # "SRR17960484_filtered_seu_2p.rds",     "g1_3",    "g1_1", "2p",
      # "SRR17960484_filtered_seu_2p.rds",     "g1_3",    "g1_5", "2p",
      # "SRR17960484_filtered_seu_2p.rds",     "s_6",      "s_7", "2p",

      "SRR13884247_filtered_seu.rds",   "g1_5-g1_7",    "g1_4-g1_0", "6p",

      "SRR17960484_filtered_seu_6p.rds",     "g1_0",      "g1_1", "6p",
      # "SRR17960484_filtered_seu_6p.rds",     "s_4",      "s_6", "6p",

      "seurat_2p_integrated_duo.rds",     "g1_1-g1_4",      "g1_9", "2p",
      "seurat_2p_integrated_duo.rds",     "g1_4",      "g1_9", "2p",
      "seurat_2p_integrated_duo.rds",     "g1_1",      "g1_9", "2p",
      "seurat_2p_integrated_duo.rds",     "g1_1",      "g1_4", "2p",

      "integrated_seu_6p_duo.rds",     "g1_4",      "g1_3", "6p",

      "SRR13884248_filtered_seu_6p.rds",   "g1_0",    "g1_2", "6p"
    ))
  ),

  # --- corresponding states: 2p ---

  tar_target(states_dictionary_2p,
    make_corresponding_states_dictionary(tibble::tribble(
      ~file_name,                       ~w_scna,    ~wo_scna, ~scna_of_interest,
      "seurat_2p_integrated_duo.rds",   "g1_1-g1_4", "g1_9",  "2p",
      "seurat_2p_integrated_duo.rds",   "g1_4",      "g1_9",  "2p",
      "seurat_2p_integrated_duo.rds",   "g1_1",      "g1_9",  "2p",
      "seurat_2p_integrated_duo.rds",   "g1_1",      "g1_4",  "2p"
    ))
  ),

  tar_target(corresponding_clusters_diffex_2p,
    find_diffex_clusters_between_corresponding_states(
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds",
      states_dictionary_2p, large_clone_comparisons,
      numbat_rds_files = numbat_rds_files, location = "all"
    )
  ),

  tar_target(corresponding_clusters_volcanos_2p,
    plot_corresponding_clusters_diffex_volcanos(
      corresponding_clusters_diffex_2p,
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"
    )
  ),

  tar_target(corresponding_clusters_enrichments_2p,
    plot_corresponding_enrichment(
      corresponding_clusters_diffex_2p,
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"
    )
  ),

  # --- corresponding states: 6p ---

  tar_target(states_dictionary_6p,
    make_corresponding_states_dictionary(tibble::tribble(
      ~file_name,                        ~w_scna,     ~wo_scna,    ~scna_of_interest,
      "SRR13884247_filtered_seu.rds",    "g1_0-g1_1", "g1_4-g1_5", "6p",
      "SRR17960484_filtered_seu_6p.rds", "g1_0",      "g1_1",      "6p"
    ))
  ),

  tar_target(corresponding_state_6p_seus,
    list(
      "output/seurat/SRR13884247_filtered_seu.rds",
      "output/seurat/SRR17960484_filtered_seu_6p.rds"
    )
  ),

  tar_target(corresponding_clusters_diffex_6p,
    find_diffex_clusters_between_corresponding_states(
      unlist(corresponding_state_6p_seus), states_dictionary_6p,
      large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"
    ),
    pattern = map(corresponding_state_6p_seus, states_dictionary_6p),
    iteration = "list"
  ),

  tar_target(corresponding_clusters_volcanos_6p,
    plot_corresponding_clusters_diffex_volcanos(
      corresponding_clusters_diffex_6p, unlist(corresponding_state_6p_seus)
    ),
    pattern = map(corresponding_clusters_diffex_6p, corresponding_state_6p_seus),
    iteration = "list"
  ),

  tar_target(corresponding_clusters_enrichments_6p,
    plot_corresponding_enrichment(
      corresponding_clusters_diffex_6p, unlist(corresponding_state_6p_seus)
    ),
    pattern = map(corresponding_clusters_diffex_6p, corresponding_state_6p_seus),
    iteration = "list"
  ),

  # --- corresponding states: combined (2p + 6p) ---

  tar_target(corresponding_clusters_diffex,
    find_diffex_clusters_between_corresponding_states(
      unlist(corresponding_seus), corresponding_states_dictionary,
      large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"
    ),
    pattern = map(corresponding_seus, corresponding_states_dictionary),
    iteration = "list"
  ),

  tar_target(corresponding_clusters_volcanos,
    plot_corresponding_clusters_diffex_volcanos(corresponding_clusters_diffex, unlist(corresponding_seus)),
    pattern = map(corresponding_clusters_diffex, corresponding_seus),
    iteration = "list"
  ),

  tar_target(corresponding_clusters_heatmaps,
    plot_corresponding_clusters_diffex_heatmaps(
      corresponding_clusters_diffex, unlist(corresponding_seus),
      corresponding_states_dictionary, large_clone_comparisons,
      numbat_rds_files = numbat_rds_files, location = "all"
    ),
    pattern = map(corresponding_clusters_diffex, corresponding_seus, corresponding_states_dictionary),
    iteration = "list"
  ),

  tar_target(corresponding_clusters_enrichments,
    plot_corresponding_enrichment(corresponding_clusters_diffex, unlist(corresponding_seus)),
    pattern = map(corresponding_clusters_diffex, corresponding_seus),
    iteration = "list"
  ),

  # --- clone cell-cycle plots ---

  tar_target(clone_cc_plots_by_scna_1q,
    clone_cc_plots_by_scna(debranched_seus_1q, scna_of_interest = "1q", large_clone_comparisons = large_clone_comparisons)
  ),

  tar_target(clone_cc_plots_by_scna_16q,
    clone_cc_plots_by_scna(debranched_seus_16q, scna_of_interest = "16q", large_clone_comparisons = large_clone_comparisons)
  ),

  # --- phase-stratified diffex (2p g1, 6p g1) ---

  tar_target(diffex_2p_g1,
    find_diffex_clones_in_phase(corresponding_seus_2p, phase = "g1", scna_of_interest = "2p", numbat_rds_files, large_clone_comparisons, location = "all"),
    pattern = map(corresponding_seus_2p),
    iteration = "list"
  ),

  tar_target(enrichment_2p_g1,
    compare_enrichment(diffex_2p_g1)
  ),

  tar_target(diffex_6p_g1,
    find_diffex_clones_in_phase(debranched_seus_6p, phase = "g1", scna_of_interest = "6p", numbat_rds_files, large_clone_comparisons, location = "all"),
    pattern = map(debranched_seus_6p),
    iteration = "list"
  ),

  tar_target(enrichment_6p_g1,
    compare_enrichment(diffex_6p_g1)
  ),

  # --- oncoprint targets ---

  tar_target(
    unfiltered_oncoprint_input_by_scna,
    make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, rb_scna_samples, n_slice = 20)
  ),

  tar_target(
    oncoprint_input_by_scna,
    filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna, oncoprint_settings)
  ),

  tar_target(
    unfiltered_oncoprint_input_by_scna_for_each_cluster,
    make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, all_diffex_clones_for_each_cluster, large_clone_comparisons, rb_scna_samples, by_cluster = TRUE, n_slice = 20)
  ),

  tar_target(
    oncoprint_input_by_scna_for_each_cluster,
    filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna_for_each_cluster, oncoprint_settings_per_cluster)
  ),

  tar_target(oncoprint_input_by_region,
    inspect_oncoprints(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans, oncoprint_input_by_scna$all)
  ),

  tar_target(oncoprint_plots,
    make_oncoprint_plots(oncoprint_input_by_scna, debranched_clone_trees, oncoprint_settings, label = "_by_clone")
  ),

  tar_target(
    oncoprint_plots_by_cluster,
    make_oncoprint_plots(oncoprint_input_by_scna_for_each_cluster, debranched_clone_trees, oncoprint_settings_per_cluster, label = "_by_cluster", p_val_threshold = 1)
  ),

  # --- oncoprint enrichment ---

  tar_target(
    oncoprint_enrich_clones_gobp,
    enrich_oncoprints(large_filter_expressions,
      cluster_dictionary, debranched_ids,
      cis_diffex_clones, trans_diffex_clones, all_diffex_clones,
      large_clone_comparisons, gene_set = "gobp"
    )
  ),

  tar_target(
    oncoprint_enrich_clones_hallmark,
    enrich_oncoprints(large_filter_expressions,
      cluster_dictionary, debranched_ids,
      cis_diffex_clones, trans_diffex_clones, all_diffex_clones,
      large_clone_comparisons, gene_set = "hallmark"
    )
  ),

  tar_target(
    oncoprint_enrich_clusters_hallmark,
    enrich_oncoprints_clusters(large_filter_expressions,
      cluster_dictionary, debranched_ids,
      cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster,
      large_clone_comparisons, gene_set = "hallmark"
    )
  ),

  tar_target(
    oncoprint_enrich_clones_plots_gobp,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_gobp,
      cis_plot_file  = "results/enrichment/cis_enrichment_plots_by_clone_gobp.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_gobp.pdf",
      cis_table_file  = "results/enrichment/cis_enrichment_tables_by_clone_gobp.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_gobp.xlsx"
    )
  ),

  tar_target(
    oncoprint_enrich_clones_plots_hallmark,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_hallmark,
      cis_plot_file  = "results/enrichment/cis_enrichment_plots_by_clone_hallmark.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_hallmark.pdf",
      cis_table_file  = "results/enrichment/cis_enrichment_tables_by_clone_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_hallmark.xlsx"
    )
  ),

  tar_target(
    oncoprint_enrich_clusters_plots_hallmark,
    compile_cis_trans_enrichment_recurrence_by_cluster(oncoprint_enrich_clusters_hallmark,
      cis_plot_file  = "results/enrichment/cis_enrichment_plots_by_cluster_hallmark.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_hallmark.pdf",
      cis_table_file  = "results/enrichment/cis_enrichment_tables_by_cluster_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_hallmark.xlsx",
      by_cluster = TRUE
    )
  ),

  # --- rod/celltype enrichment ---

  tar_target(rod_rich_samples,
    score_samples_for_rod_enrichment(numbat_rds_files),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),

  tar_target(table_06,
    rod_rich_samples
  ),

  tar_target(celltype_rich_samples,
    score_samples_for_celltype_enrichment(unfiltered_seus, final_seus, celltype_markers),
    pattern = map(unfiltered_seus),
    iteration = "list"
  ),

  tar_target(stachelek_score_plots,
    score_stachelek(final_seus, oncoprint_input_by_scna),
    pattern = map(final_seus),
    iteration = "list"
  ),

  # --- pseudobulks ---

  tar_target(whole_pseudobulks,
    score_whole_pseudobulks(numbat_rds_files, subtype_markers),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),

  tar_target(
    unfiltered_derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(unfiltered_seus)
  ),

  tar_target(
    filtered_derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(final_seus)
  ),

  # --- figures 09 / 10 ---

  tar_target(fig_09,
    plot_fig_09_10(
      corresponding_seus_2p, corresponding_seus,
      corresponding_clusters_diffex, corresponding_clusters_enrichments,
      recurrence_threshold = 3, plot_path = "results/fig_09.pdf",
      widths = rep(4, 3), heights = rep(8, 3),
      common_seus = c("SRR13884248_filtered_seu_2p.rds", "SRR17960484_filtered_seu_2p.rds")
    )
  ),

  tar_target(fig_10,
    plot_fig_09_10(
      corresponding_seus_6p, corresponding_seus,
      corresponding_clusters_diffex, corresponding_clusters_enrichments,
      recurrence_threshold = 2, plot_path = "results/fig_10.pdf",
      widths = rep(4, 3), heights = c(12, 4, 12),
      common_seus = c("SRR13884247_filtered_seu_6p.rds", "SRR17960484_filtered_seu_6p.rds")
    )
  )

)
