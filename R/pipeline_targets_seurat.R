# Seurat processing targets: filtering, debranching, clustering, heatmap collages,
# regression diagnostics, hypoxia, integrated seus, and metadata DB.
# Defines: pipeline_targets_seurat (spliced into tar_plan in _targets.R)

pipeline_targets_seurat <- c(

  # --- file-tracking: debranched SCNA seu paths ---

  # format = "file" so that changes to the RDS files on disk invalidate
  # downstream targets (e.g. factor heatmaps, diffex, clone pearls).
  tarchetypes::tar_files(debranched_seus_1q,
    c(
      "SRR13884246" = "output/seurat/SRR13884246_branch_5_filtered_seu_1q.rds",
      "SRR13884249" = "output/seurat/SRR13884249_filtered_seu_1q.rds",
      "SRR14800534" = "output/seurat/SRR14800534_filtered_seu_1q.rds",
      "SRR14800535" = "output/seurat/SRR14800535_filtered_seu_1q.rds",
      "SRR14800536" = "output/seurat/SRR14800536_filtered_seu_1q.rds"
    )
  ),

  tarchetypes::tar_files(debranched_seus_2p,
    c(
      "SRR13884246" = "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds",
      "SRR13884247" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
      "SRR13884248" = "output/seurat/SRR13884248_filtered_seu_2p.rds",
      # "SRR13884249" = "output/seurat/SRR13884249_filtered_seu_2p.rds",
      # "SRR17960481" = "output/seurat/SRR17960481_filtered_seu.rds",
      "SRR17960484" = "output/seurat/SRR17960484_filtered_seu_2p.rds" # too many changes at once to conclude anything
    )
  ),

  tarchetypes::tar_files(debranched_seus_6p,
    c(
      "SRR13884247" = "output/seurat/SRR13884247_filtered_seu.rds",
      "SRR13884248" = "output/seurat/SRR13884248_filtered_seu_6p.rds",
      "SRR17960484" = "output/seurat/SRR17960484_filtered_seu_6p.rds",
      "SRR27187899" = "output/seurat/SRR27187899_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(debranched_seus_16q,
    c(
      "SRR14800534" = "output/seurat/SRR14800534_filtered_seu_16q.rds",
      "SRR14800535" = "output/seurat/SRR14800535_filtered_seu_16q.rds",
      "SRR14800536" = "output/seurat/SRR14800536_filtered_seu_16q.rds"
    )
  ),

  # --- file-tracking: integrated seu paths ---

  tarchetypes::tar_files(integrated_seus_1q,
    c(
      "SRR13884249_filtered_seu.rds" = "output/seurat/integrated_1q/SRR13884249_integrated_1q_filtered_seu.rds",
      "SRR14800534_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800534_integrated_1q_filtered_seu.rds",
      "SRR14800535_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800535_integrated_1q_filtered_seu.rds",
      "SRR14800536_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800536_integrated_1q_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(integrated_seus_16q,
    c(
      "SRR14800534_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800534_integrated_16q_filtered_seu.rds",
      "SRR14800535_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800535_integrated_16q_filtered_seu.rds",
      "SRR14800536_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800536_integrated_16q_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(integrated_seus_2p,
    c(
      "SRR13884248_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRR13884248_integrated_2p_filtered_seu.rds",
      "SRR17960484_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRR17960484_integrated_2p_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(integrated_seus_6p,
    c(
      "SRR13884248_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRR13884248_integrated_6p_filtered_seu.rds",
      "SRR17960484_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRR17960484_integrated_6p_filtered_seu.rds"
    )
  ),

  list(

    # --- per-sample Seurat processing ---

    tar_target(filtered_large_plot_files,
      make_numbat_plot_files(seus, numbat_rds_files, cluster_dictionary, large_filter_expressions, large_clone_simplifications, extension = "_filtered"),
      pattern = map(seus, numbat_rds_files),
      iteration = "list"
    ),

    tar_target(unfiltered_seus,
      prep_unfiltered_seu(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, extension = "_unfiltered"),
      pattern = map(numbat_rds_files),
      iteration = "list"
    ),

    tar_target(filtered_seus,
      filter_cluster_save_seu(numbat_rds_files, seus, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove, extension = "", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
      pattern = map(numbat_rds_files),
      iteration = "list"
    ),

    tar_target(
      final_seus,
      set_final_seus(interesting_samples)
    ),

    tar_target(
      cc_plots_wo_arms,
      plot_phase_wo_arm(final_seus)
    ),

    # --- debranched seus ---

    tar_file(debranched_seu_files,
      c(
        "SRR13884242"          = "output/seurat/SRR13884242_filtered_seu.rds",
        "SRR13884243"          = "output/seurat/SRR13884243_filtered_seu.rds",
        "SRR13884246_branch_5" = "output/seurat/SRR13884246_branch_5_filtered_seu.rds",
        "SRR13884246_branch_6" = "output/seurat/SRR13884246_branch_6_filtered_seu.rds",
        "SRR13884247_branch_6" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
        "SRR13884247_branch_4" = "output/seurat/SRR13884247_branch_4_filtered_seu.rds",
        "SRR13884247_branch_5" = "output/seurat/SRR13884247_branch_5_filtered_seu.rds",
        "SRR13884248"          = "output/seurat/SRR13884248_filtered_seu.rds",
        "SRR13884249"          = "output/seurat/SRR13884249_filtered_seu.rds",
        "SRR14800534"          = "output/seurat/SRR14800534_filtered_seu.rds",
        "SRR14800535"          = "output/seurat/SRR14800535_filtered_seu.rds",
        "SRR14800536"          = "output/seurat/SRR14800536_filtered_seu.rds",
        "SRR14800540_branch_2" = "output/seurat/SRR14800540_branch_2_filtered_seu.rds",
        "SRR14800540_branch_3" = "output/seurat/SRR14800540_branch_3_filtered_seu.rds",
        "SRR14800541_branch_4" = "output/seurat/SRR14800541_branch_4_filtered_seu.rds",
        "SRR14800541_branch_7" = "output/seurat/SRR14800541_branch_7_filtered_seu.rds",
        "SRR14800543_branch_3" = "output/seurat/SRR14800543_branch_3_filtered_seu.rds",
        "SRR14800543_branch_4" = "output/seurat/SRR14800543_branch_4_filtered_seu.rds",
        "SRR17960481"          = "output/seurat/SRR17960481_filtered_seu.rds",
        "SRR17960484"          = "output/seurat/SRR17960484_filtered_seu.rds",
        "SRR27187899"          = "output/seurat/SRR27187899_filtered_seu.rds",
        "SRR27187902_branch_3" = "output/seurat/SRR27187902_branch_3_filtered_seu.rds",
        "SRR27187902_branch_4" = "output/seurat/SRR27187902_branch_4_filtered_seu.rds"
      )
    ),

    tar_target(debranched_seus,
      set_names(debranched_seu_files, str_remove(fs::path_file(debranched_seu_files), "_filtered_seu.rds"))
    ),

    tar_target(overall_seus,
      merge_orders(debranched_seus, final_seus)
    ),

    # --- scna-stratified seu collections ---

    tar_target(scna_seus,
      unlist(list(
        "1q"  = debranched_seus_1q,
        "2p"  = debranched_seus_2p,
        "6p"  = debranched_seus_6p,
        "16q" = debranched_seus_16q
      ))
    ),

    # Populate SQLite metadata tables for all SCNA-stratified Seurat objects.
    # Upserts four tables: seurat_objects, cell_metadata, cluster_composition,
    # qc_metrics. Re-runs automatically when any debranched_seus_* file changes.
    tar_target(seu_metadata_db,
      bulk_extract_seu_metadata(scna_seus, sqlite_path = "batch_hashes.sqlite"),
      deployment = "main"  # DB writes must run on the main process, not a worker
    ),

    tar_target(resolution_dictionary,
      read_resolution_dictionary(sqlite_path = "batch_hashes.sqlite"),
      deployment = "main"
    ),

    tar_target(chosen_resolution_seus,
      assign_designated_phase_clusters(scna_seus, cluster_orders, resolution_dictionary),
      pattern = map(scna_seus),
      iteration = "list"
    ),

    # --- clustrees ---

    # --- recurrent heatmaps ---

    tar_target(
      filtered_recurrent_genes,
      pull_common_markers(final_seus, mps[["Cancer"]])
    ),

    tar_target(filtered_recurrent_heatmaps,
      heatmap_marker_genes(
        final_seus,
        filtered_recurrent_genes,
        subtype_markers,
        "filtered_",
        marker_col = "SCT_snn_res.0.2",
        group.by = c("SCT_snn_res.0.2", "Phase", "scna"),
        col_arrangement = c("SCT_snn_res.0.2", "Phase", "scna")
      ),
      pattern = map(final_seus),
      iteration = "list"
    ),

    tar_target(
      filtered_recurrent_heatmap_file,
      qpdf::pdf_combine(filtered_recurrent_heatmaps, "results/filtered_recurrent_heatmaps.pdf")
    ),

    tar_target(
      combined_recurrent_filtered_heatmap,
      heatmap_marker_genes(
        "output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds",
        filtered_recurrent_genes,
        subtype_markers,
        "filtered_",
        label = "SRR14800534_SRR14800535_SRR14800536",
        marker_col = "SCT_snn_res.0.2",
        group.by = c("integrated_snn_res.0.2", "Phase", "scna"),
        col_arrangement = c("integrated_snn_res.0.2", "Phase", "scna")
      )
    ),

    tar_target(
      regressed_recurrent_genes,
      pull_common_markers(regressed_seus, mps[["Cancer"]])
    ),

    tar_target(regressed_recurrent_heatmaps,
      heatmap_marker_genes(
        regressed_seus,
        regressed_recurrent_genes,
        subtype_markers,
        "regressed_",
        marker_col = "SCT_snn_res.0.2",
        group.by = c("SCT_snn_res.0.2", "Phase", "scna"),
        col_arrangement = c("SCT_snn_res.0.2", "scna")
      ),
      pattern = map(regressed_seus),
      iteration = "list"
    ),

    tar_target(regressed_recurrent_heatmap_file,
      qpdf::pdf_combine(regressed_recurrent_heatmaps, "results/regressed_recurrent_heatmaps.pdf")
    ),

    # --- regression ---

    tar_target(regressed_seus,
      regress_filtered_seu(final_seus),
      pattern = map(final_seus),
      iteration = "list"
    ),

    tar_target(effect_of_filtering,
      plot_effect_of_filtering(unfiltered_seus, final_seus),
      pattern = map(unfiltered_seus),
      iteration = "list"
    ),

    tar_target(regression_effect_plots,
      plot_effect_of_regression(final_seus, regressed_seus, width = 18, height = 12),
      pattern = map(final_seus, regressed_seus),
      iteration = "list"
    ),

    # regression diagnostics
    tar_target(fig_03_05, {
      plot_paths <- na.omit(unlist(regression_effect_plots))
      plot_paths <- plot_paths[nzchar(plot_paths)]
      out_path <- "results/fig_03_05.pdf"
      if (length(plot_paths) == 0) {
        pdf(out_path)
        plot.new()
        text(0.5, 0.5, "No valid regression effect plots")
        dev.off()
        out_path
      } else {
        qpdf::pdf_combine(plot_paths, out_path)
      }
    }),

    tar_target(regression_ora_plots,
      ora_effect_of_regression(final_seus, regressed_seus),
      pattern = map(final_seus, regressed_seus),
      iteration = "list"
    ),

    tar_target(cin_score_plots,
      score_chrom_instability(unfiltered_seus),
      pattern = map(unfiltered_seus),
      iteration = "list"
    ),

    # --- QC / numbat heatmaps ---

    tar_target(fig_s03a_plots,
      make_numbat_heatmaps(original_seus, numbat_rds_files, p_min = 0.9, line_width = 0.1, extension = "_unfiltered"),
      pattern = map(original_seus),
      iteration = "list"
    ),

    tar_target(fig_s03a, qpdf::pdf_combine(unlist(fig_s03a_plots), "results/unfiltered_heatmaps.pdf")),

    tar_target(fig_s13,
      make_numbat_heatmaps(final_seus, numbat_rds_files, p_min = 0.5, line_width = 0.1, extension = "_filtered"),
      pattern = map(final_seus),
      iteration = "list",
      cue = tar_cue("always")  # force re-render each run: heatmap layout depends on
                                # session graphics device state, not captured by hash
    ),

    tar_target(filtered_numbat_heatmaps_file, qpdf::pdf_combine(map_chr(fig_s13, 1), "results/filtered_heatmaps.pdf")),

    tar_target(filtered_large_scna_prob_file, {
      scna_var_paths <- na.omit(map_chr(fig_s13, 2))
      qpdf::pdf_combine(scna_var_paths, "results/filtered_scna_probabilities.pdf")
    }),

    # --- heatmap collages ---

    tar_target(heatmap_collages,
      plot_seu_marker_heatmap_all_resolutions(debranched_seus, numbat_rds_files, large_clone_simplifications),
      pattern = map(debranched_seus),
      iteration = "list"
    ),

    tar_target(heatmap_collages_6p,
      plot_seu_marker_heatmap_all_resolutions(debranched_seus_6p, numbat_rds_files, large_clone_simplifications),
      pattern = map(debranched_seus_6p),
      iteration = "list"
    ),

    tar_target(annotated_heatmap_collages,
      plot_seu_marker_heatmap(debranched_seus, cluster_orders, numbat_rds_files, large_clone_simplifications),
      pattern = map(debranched_seus),
      iteration = "list"
    ),

    tar_target(
      collage_compilation,
      qpdf::pdf_combine(annotated_heatmap_collages, "results/heatmap_collages.pdf")
    ),

    tar_target(
      collage_compilation_all_resolutions,
      qpdf::pdf_combine(heatmap_collages, "results/heatmap_collages_all_resolutions.pdf")
    ),

    # --- hypoxia ---

    tar_target(hypoxia_seus,
      load_and_save_hypoxia_score(filtered_seus),
      pattern = map(filtered_seus),
      iteration = "list"
    ),

    tar_target(
      hypoxia_score_plots,
      plot_hypoxia_score(hypoxia_seus, threshold = 0.5),
      pattern = map(hypoxia_seus),
      iteration = "list"
    ),

    tar_target(seus_low_hypoxia,
      subset_seu_by_expression(hypoxia_seus, run_hypoxia_clustering = TRUE, hypoxia_expr = "hypoxia_score <= 0.5", slug = "hypoxia_low"),
      pattern = map(hypoxia_seus),
      iteration = "list"
    ),

    tar_target(seus_high_hypoxia,
      subset_seu_by_expression(hypoxia_seus, run_hypoxia_clustering = TRUE, hypoxia_expr = "hypoxia_score > 0.5", slug = "hypoxia_high"),
      pattern = map(hypoxia_seus),
      iteration = "list"
    ),

    tar_target(heatmap_collages_hypoxia,
      plot_seu_marker_heatmap(hypoxia_seus, nb_paths = numbat_rds_files, clone_simplifications = large_clone_simplifications, tmp_plot_path = TRUE),
      pattern = map(hypoxia_seus),
      iteration = "list"
    ),

    # TODO: heatmap_collages_hypoxia_low and heatmap_collages_hypoxia_high need to be
    # updated to use seus_low_hypoxia / seus_high_hypoxia respectively.

    tar_target(hypoxia_effect_plots,
      list(
        qpdf::pdf_combine(unlist(hypoxia_score_plots), "01_hypoxia_score_plots.pdf"),
        qpdf::pdf_combine(unlist(heatmap_collages_hypoxia), "02_hypoxia_heatmap.pdf")
        # TODO: add low/high hypoxia heatmaps once those targets are fixed
      ),
    ),

    tar_target(silhouette_plots,
      calc_silhouette(debranched_seus),
      pattern = map(debranched_seus),
      iteration = "list"
    ),

    # --- integrated seus ---

    tar_target(
      integrated_seu_16q,
      readRDS("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds")
    ),

    tar_target(integrated_seus,
      list(
        "diploid_v_16q"    = "output/seurat/SRR1_diploid_v_16q_filtered_seu.rds",
        "16q_v_16q-1q"     = "output/seurat/SRR1_16q_v_16q-1q_filtered_seu.rds",
        "diploid_v_16q-1q" = "output/seurat/SRR1_diploid_v_16q-1q_filtered_seu.rds"
      )
    ),

    tar_target(annotated_integrated_heatmap_collages,
      plot_seu_marker_heatmap_integrated(integrated_seus),
      pattern = map(integrated_seus),
      iteration = "list"
    )

  ),

  # clustree_{id} — one per debranched sample/branch (static branching so
  # targets are named e.g. clustree_SRR13884242, clustree_SRR13884246_branch_5)
  # seu_file_{id} provides per-sample file tracking so only the affected
  # clustree reruns when a single Seurat file changes on disk.
  {
    .clustree_map <- tarchetypes::tar_map(
      unlist = FALSE,
      values = debranched_map_values,
      names  = "id",
      tarchetypes::tar_file(seu_file, seu_path),
      tar_target(clustree,
        make_clustrees_for_sample(seu_file, mylabel = id, assay = "SCT")
      )
    )
    .clustree_combined <- tarchetypes::tar_combine(
      clustrees,
      .clustree_map[["clustree"]],
      command = list(!!!.x)
    )
    list(
      .clustree_map,
      .clustree_combined,
      tarchetypes::tar_files(
        clustree_compilation,
        qpdf::pdf_combine(unlist(clustrees), "results/clustrees.pdf"),
        format = "file"
      )
    )
  },

  # hypoxia_seus_1q/2p/6p/16q — subset seus_low_hypoxia to each SCNA's samples
  tarchetypes::tar_map(
    values = scna_map_values[, "scna", drop = FALSE],
    names  = "scna",
    tar_target(hypoxia_seus,
      str_subset(unlist(seus_low_hypoxia), str_c(rb_scna_samples[[scna]], collapse = "|"))
    )
  ),

  # integrated_seu_low_hypoxia_1q/2p/6p/16q — integration within each SCNA at low hypoxia
  tarchetypes::tar_map(
    values = tibble::tibble(
      scna        = c("1q", "2p", "6p", "16q"),
      hypoxia_sym = rlang::syms(c("hypoxia_seus_1q", "hypoxia_seus_2p",
                                  "hypoxia_seus_6p", "hypoxia_seus_16q"))
    ),
    names = "scna",
    tar_target(integrated_seu_low_hypoxia,
      integration_by_scna_clones(
        hypoxia_sym, scna_of_interest = scna,
        large_clone_comparisons, filter_expr = "hypoxia_score <= 0.5"
      )
    )
  ),

  list(

    tar_target(hypoxia_low_integrated_heatmap_collages,
      plot_seu_marker_heatmap_integrated(integrated_seu_low_hypoxia_1q)
    ),

    tar_target(hypoxia_low_integrated_heatmap_collages0,
      plot_integrated_1q_fig_low_hypoxia(integrated_seu_low_hypoxia_1q)
    )

  )

)
