# Seurat processing targets: filtering, debranching, clustering, heatmap collages,
# regression diagnostics, hypoxia, integrated seus, and metadata DB.
# Defines: pipeline_targets_seurat (spliced into tar_plan in _targets.R)

# Crew controller resource assignments for targets
.light_resources <- tar_resources(crew = tar_resources_crew(controller = "light"))
.heavy_resources <- tar_resources(crew = tar_resources_crew(controller = "heavy"))

pipeline_targets_seurat <- c(

  # --- file-tracking: debranched SCNA seu paths ---

  # format = "file" so that changes to the RDS files on disk invalidate
  # downstream targets (e.g. factor heatmaps, diffex, clone pearls).
  tarchetypes::tar_files(debranched_seus_1q,
    c(
      "SRX10264523" = "output/seurat/SRX10264523_branch_5_filtered_seu_1q.rds",
      "SRX10264526" = "output/seurat/SRX10264526_filtered_seu_1q.rds",
      "SRX11133594" = "output/seurat/SRX11133594_filtered_seu_1q.rds",
      "SRX11133593" = "output/seurat/SRX11133593_filtered_seu_1q.rds",
      "SRX11133592" = "output/seurat/SRX11133592_filtered_seu_1q.rds",
      # "SRX10831287" = "output/seurat/SRX10831287_filtered_seu_1q.rds"
    )
  ),

  tarchetypes::tar_files(debranched_seus_2p,
    c(
      "SRX10264523" = "output/seurat/SRX10264523_branch_5_filtered_seu_2p.rds",
      "SRX10264524" = "output/seurat/SRX10264524_branch_6_filtered_seu.rds",
      "SRX10264525" = "output/seurat/SRX10264525_filtered_seu_2p.rds",
      # "SRX10264526" = "output/seurat/SRX10264526_filtered_seu_2p.rds",
      # "SRX14116947" = "output/seurat/SRX14116947_filtered_seu.rds",
      "SRX14116944" = "output/seurat/SRX14116944_filtered_seu_2p.rds" # too many changes at once to conclude anything
    )
  ),

  tarchetypes::tar_files(debranched_seus_6p,
    c(
      "SRX10264524" = "output/seurat/SRX10264524_filtered_seu.rds",
      "SRX10264525" = "output/seurat/SRX10264525_filtered_seu_6p.rds",
      "SRX14116944" = "output/seurat/SRX14116944_filtered_seu_6p.rds",
      "SRX22868105" = "output/seurat/SRX22868105_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(debranched_seus_16q,
    c(
      "SRX11133594" = "output/seurat/SRX11133594_filtered_seu_16q.rds",
      "SRX11133593" = "output/seurat/SRX11133593_filtered_seu_16q.rds",
      "SRX11133592" = "output/seurat/SRX11133592_filtered_seu_16q.rds"
    )
  ),

  # --- file-tracking: integrated seu paths ---

  tarchetypes::tar_files(integrated_seus_1q,
    c(
      "SRX10264526_filtered_seu.rds" = "output/seurat/integrated_1q/SRX10264526_integrated_1q_filtered_seu.rds",
      "SRX11133594_filtered_seu.rds" = "output/seurat/integrated_1q/SRX11133594_integrated_1q_filtered_seu.rds",
      "SRX11133593_filtered_seu.rds" = "output/seurat/integrated_1q/SRX11133593_integrated_1q_filtered_seu.rds",
      "SRX11133592_filtered_seu.rds" = "output/seurat/integrated_1q/SRX11133592_integrated_1q_filtered_seu.rds"
      # SRX10831287 excluded: integrated seu not yet produced
    )
  ),

  tarchetypes::tar_files(integrated_seus_16q,
    c(
      "SRX11133594_filtered_seu.rds" = "output/seurat/integrated_16q/SRX11133594_integrated_16q_filtered_seu.rds",
      "SRX11133593_filtered_seu.rds" = "output/seurat/integrated_16q/SRX11133593_integrated_16q_filtered_seu.rds",
      "SRX11133592_filtered_seu.rds" = "output/seurat/integrated_16q/SRX11133592_integrated_16q_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(integrated_seus_2p,
    c(
      "SRX10264525_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRX10264525_integrated_2p_filtered_seu.rds",
      "SRX14116944_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRX14116944_integrated_2p_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_files(integrated_seus_6p,
    c(
      "SRX10264525_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRX10264525_integrated_6p_filtered_seu.rds",
      "SRX14116944_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRX14116944_integrated_6p_filtered_seu.rds"
    )
  ),

  list(

    # --- per-sample Seurat processing ---

    tar_target(filtered_large_plot_files,
      # Build filtered large-format Numbat plot bundles per sample.
      make_numbat_plot_files(seus_interesting, numbat_rds_files, cluster_dictionary, large_filter_expressions, large_clone_simplifications, extension = "_filtered"),
      pattern = map(seus_interesting, numbat_rds_files),
      iteration = "list"
    ),

    tar_target(unfiltered_seus,
      # First stage where scna metadata is injected into Seurat objects.
      prep_unfiltered_seu(numbat_rds_files, cluster_dictionary_per_sample, large_clone_simplifications_per_sample, large_filter_expressions_per_sample, extension = "_unfiltered"),
      pattern = map(numbat_rds_files, cluster_dictionary_per_sample, large_clone_simplifications_per_sample, large_filter_expressions_per_sample),
      iteration = "list",
      error = "null",
      cue = tar_cue(command = FALSE, depend = FALSE)
    ),

    tar_target(filter_inspection_metadata,
      # Lightweight per-cell metadata tibble for threshold sweeping.
      # Reads unfiltered_seus (no SCTransform/UMAP); invalidated only when Seurat content changes.
      extract_filter_metadata(
        unfiltered_seus,
        cluster_dictionary_per_sample,
        cells_to_remove,
        large_clone_simplifications_per_sample
      ),
      pattern = map(unfiltered_seus, cluster_dictionary_per_sample, large_clone_simplifications_per_sample),
      iteration = "list"
    ),

    tar_target(filtered_seus,
      # Filtered counterpart also computes/adds scna metadata at cell level.
      filter_cluster_save_seu(
        numbat_rds_files, unfiltered_seus,
        cluster_dictionary_per_sample, large_clone_simplifications_per_sample,
        filter_expressions = NULL, cells_to_remove,
        extension = "",
        leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"
      ),
      pattern = map(numbat_rds_files, cluster_dictionary_per_sample, large_clone_simplifications_per_sample),
      iteration = "list",
      resources = .heavy_resources,
      error = "null",
      cue = tar_cue(command = FALSE, depend = FALSE)  # re-run when unfiltered_seus changes, but not when filter_inspection_metadata changes
    ),

    tar_target(filtered_seus_with_phase,
      # Calculate Phase (cell cycle) for each filtered_seu before integration
      add_phase_to_filtered_seu(filtered_seus),
      pattern = map(filtered_seus),
      iteration = "list",
      error = "null",
      cue = tar_cue(command = FALSE)
    ),

    tar_target(rod_low_sample_ids,
      dplyr::bind_rows(table_rod_rich_samples) |>
        dplyr::filter(rod_rich == 0) |>
        dplyr::pull(sample_id)
    ),

    tar_target(ks_diploid_seu,
      assemble_diploid_seu(
        seus_low_hypoxia[grepl(
          paste(rod_low_sample_ids, collapse = "|"),
          unlist(seus_low_hypoxia)
        ) & !grepl(
          "SRX10031194|SRX10264517|SRX10264518|SRX10264523|SRX14116946",
          unlist(seus_low_hypoxia)
        )],
        integrate = TRUE
      ),
      format = "file"
    ),

    tar_target(diploid_seu,
      "output/seurat/diploid_subsets/diploid_seu.rds",
      format = "file"
    ),

    tar_target(diploid_seu_umap_plots,
      plot_diploid_seu_umaps(diploid_seu, celltype_markers),
      format = "file"
    ),

    tar_target(filtered_seus_nb_filtered,
      # Seurat objects annotated with clone/SCNA metadata from numbat_sridhar_filtered RDS files.
      filter_cluster_save_seu(
        numbat_rds_filtered_files, unfiltered_seus,
        cluster_dictionary, large_clone_simplifications,
        filter_expressions = NULL, cells_to_remove,
        extension = "_nb_filtered",
        leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"
      ),
      pattern = map(numbat_rds_filtered_files),
      iteration = "list",
      resources = .heavy_resources,
      error = "null",
      cue = tar_cue(command = FALSE, depend = FALSE)
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
        "SRX10264519"          = "output/seurat/SRX10264519_filtered_seu.rds",
        "SRX10264520"          = "output/seurat/SRX10264520_filtered_seu.rds",
        "SRX10264523_branch_5" = "output/seurat/SRX10264523_branch_5_filtered_seu.rds",
        "SRX10264523_branch_6" = "output/seurat/SRX10264523_branch_6_filtered_seu.rds",
        "SRX10264524_branch_6" = "output/seurat/SRX10264524_branch_6_filtered_seu.rds",
        "SRX10264524_branch_4" = "output/seurat/SRX10264524_branch_4_filtered_seu.rds",
        "SRX10264524_branch_5" = "output/seurat/SRX10264524_branch_5_filtered_seu.rds",
        "SRX10264525"          = "output/seurat/SRX10264525_filtered_seu.rds",
        "SRX10264526"          = "output/seurat/SRX10264526_filtered_seu.rds",
        "SRX11133594"          = "output/seurat/SRX11133594_filtered_seu.rds",
        "SRX11133593"          = "output/seurat/SRX11133593_filtered_seu.rds",
        "SRX11133592"          = "output/seurat/SRX11133592_filtered_seu.rds",
        "SRX11133588_branch_2" = "output/seurat/SRX11133588_branch_2_filtered_seu.rds",
        "SRX11133588_branch_3" = "output/seurat/SRX11133588_branch_3_filtered_seu.rds",
        "SRX11133587_branch_4" = "output/seurat/SRX11133587_branch_4_filtered_seu.rds",
        "SRX11133587_branch_7" = "output/seurat/SRX11133587_branch_7_filtered_seu.rds",
        "SRX11133585_branch_3" = "output/seurat/SRX11133585_branch_3_filtered_seu.rds",
        "SRX11133585_branch_4" = "output/seurat/SRX11133585_branch_4_filtered_seu.rds",
        "SRX14116947"          = "output/seurat/SRX14116947_filtered_seu.rds",
        "SRX14116944"          = "output/seurat/SRX14116944_filtered_seu.rds",
        "SRX22868105"          = "output/seurat/SRX22868105_filtered_seu.rds",
        "SRX22868102_branch_3" = "output/seurat/SRX22868102_branch_3_filtered_seu.rds",
        "SRX22868102_branch_4" = "output/seurat/SRX22868102_branch_4_filtered_seu.rds"
      )
    ),

    tar_target(debranched_seus,
      set_names(debranched_seu_files, str_remove(fs::path_file(debranched_seu_files), "_filtered_seu.rds"))
    ),

    tar_target(overall_seus,
      merge_orders(seus_low_hypoxia, final_seus)
    ),

    # --- scna-stratified seu collections ---

    tar_target(scna_seus,
      # Canonical union of SCNA-stratified Seurat paths used by DB extraction.
      unlist(list(
        "1q"  = hypoxia_seus_1q,
        "2p"  = hypoxia_seus_2p,
        "6p"  = hypoxia_seus_6p,
        "16q" = hypoxia_seus_16q
      ))
    ),

    # Collect all available Seurat RDS files (unfiltered + filtered from all samples)
    tar_target(all_seu_files,
      setNames(as.character(fs::dir_ls("output/seurat/", regexp = "/SRX[0-9]+_seu.rds")),
               stringr::str_extract(fs::dir_ls("output/seurat/", regexp = "/SRX[0-9]+_seu.rds"), "SRX[0-9]+")),
      deployment = "main"
    ),

    # Populate SQLite metadata tables for all Seurat objects.
    # Upserts eight tables: seurat_objects, cell_metadata, cell_qc_values,
    # cluster_composition, cluster_markers, qc_metrics, hashes, cluster_orders.
    # Re-runs automatically when any .rds file in output/seurat/ changes.
    tar_target(seu_metadata_db,
      bulk_extract_seu_metadata(all_seu_files, sqlite_path = "batch_hashes.sqlite"),
      deployment = "main"  # DB writes must run on the main process, not a worker
    ),

    # Derive _hypoxia_low_seu.rds aliases from existing _filtered_seu_*.rds entries
    # so that assign_designated_phase_clusters finds resolutions for hypoxia paths.
    tar_target(hypoxia_resolution_aliases, {
      existing <- read_resolution_dictionary(sqlite_path = "batch_hashes.sqlite")
      if (length(existing) == 0) return(invisible(0L))
      tumor_ids <- stringr::str_extract(names(existing), "SR[RX][0-9]+")
      hypoxia_keys <- paste0(tumor_ids, "_hypoxia_low_seu.rds")
      alias_df <- data.frame(
        file_id     = hypoxia_keys,
        resolution  = unlist(existing),
        stringsAsFactors = FALSE
      )
      alias_df <- alias_df[!alias_df$file_id %in% names(existing), , drop = FALSE]
      if (nrow(alias_df) > 0) upsert_resolution_dictionary(alias_df)
      invisible(nrow(alias_df))
    }, deployment = "main", cue = tar_cue("always")),

    tar_target(resolution_dictionary, {
      hypoxia_resolution_aliases
      read_resolution_dictionary(sqlite_path = "batch_hashes.sqlite")
    }, deployment = "main"),

    tar_target(chosen_resolution_seus,
      # Attach designated phase cluster resolutions pulled from sqlite metadata.
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
        "output/seurat/SRX11133594_SRX11133593_SRX11133592_seu.rds",
        filtered_recurrent_genes,
        subtype_markers,
        "filtered_",
        label = "SRX11133594_SRX11133593_SRX11133592",
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
      iteration = "list",
      error = "null",
      cue = tar_cue(command = FALSE)
    ),

    tar_target(effect_of_filtering,
      {
        path <- unlist(unfiltered_seus)
        if (is.na(path)) return(NULL)
        sample_id <- stringr::str_extract(path, "SRX[0-9]+")
        filtered_path <- unlist(filtered_seus)[grepl(sample_id, unlist(filtered_seus))]
        if (length(filtered_path) == 0) filtered_path <- NULL
        lh_paths <- unlist(seus_low_hypoxia)
        low_hypoxia_path <- lh_paths[grepl(sample_id, lh_paths)]
        if (length(low_hypoxia_path) == 0) low_hypoxia_path <- NULL
        hh_paths <- unlist(seus_high_hypoxia)
        high_hypoxia_path <- hh_paths[grepl(sample_id, hh_paths)]
        if (length(high_hypoxia_path) == 0) high_hypoxia_path <- NULL
        plot_effect_of_filtering(path, filtered_path, cluster_dictionary = cluster_dictionary,
                                 low_hypoxia_seu_path = low_hypoxia_path,
                                 high_hypoxia_seu_path = high_hypoxia_path)
      },
      pattern = map(unfiltered_seus),
      iteration = "list",
      error = "null",
      cue = tar_cue(depend = FALSE)
    ),

    tar_target(effect_of_filtering_collated,
      {
        paths <- unlist(effect_of_filtering)
        paths <- paths[!sapply(paths, is.null) & file.exists(paths)]
        out_path <- "results/effect_of_filtering_markers.pdf"
        qpdf::pdf_combine(paths, out_path)
      },
      format = "file"
    ),

    # Cell counts at each filtering stage.
    # Stages 1-2 come from filter_inspection_metadata (no Seurat reads).
    # Stages 3-4 query batch_hashes.sqlite, populated by add_hash_metadata
    # when filtered_seus / seus_low_hypoxia are computed.
    tar_target(filtering_cell_counts_table,
      generate_filtering_cell_counts(filtered_seus, seus_low_hypoxia, filter_inspection_metadata),
      format = "file"
    ),

    tar_target(regression_effect_plots,
      plot_effect_of_regression(final_seus, regressed_seus, width = 18, height = 12),
      pattern = map(final_seus, regressed_seus),
      iteration = "list"
    ),

    # regression diagnostics
    tar_target(fig_regression_diagnostics, {
      plot_paths <- na.omit(unlist(regression_effect_plots))
      plot_paths <- plot_paths[nzchar(plot_paths)]
      out_path <- "results/fig_regression_diagnostics.pdf"
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

    tar_target(numbat_heatmap_plots_unfiltered,
      make_numbat_heatmaps(
        unfiltered_seus, numbat_rds_files,
        p_min = 0.9, line_width = 0.1, extension = "_unfiltered",
        show_segment_names_on_x = TRUE, filter_midline = FALSE
      ),
      pattern = map(unfiltered_seus),
      iteration = "list"
    ),

    tar_target(numbat_heatmap_plots_subset,
      make_numbat_heatmaps(
        filtered_seus, numbat_rds_files,
        p_min = 0.9, line_width = 0.1, extension = "_subset",
        show_segment_names_on_x = TRUE, filter_midline = FALSE
      ),
      pattern = map(filtered_seus),
      iteration = "list"
    ),

    tar_target(numbat_heatmaps_unfiltered_pdf, {
      paths <- unlist(numbat_heatmap_plots_unfiltered)
      paths <- paths[!is.na(paths)]
      qpdf::pdf_combine(paths, "results/unfiltered_heatmaps.pdf")
    }),

    tar_target(numbat_heatmap_plots_low_hypoxia,
      make_numbat_heatmaps(
        seus_low_hypoxia, numbat_rds_files,
        p_min = 0.9, line_width = 0.1, extension = "_low_hypoxia",
        show_segment_names_on_x = TRUE,
        numbat_rds_filtered_files = numbat_rds_files,
        filter_midline = FALSE
      ),
      pattern = map(seus_low_hypoxia),
      iteration = "list"
    ),

    tar_target(fig_numbat_heatmaps, {
      paths <- unlist(numbat_heatmap_plots_low_hypoxia)
      paths <- paths[!is.na(paths)]
      qpdf::pdf_combine(paths, "results/low_hypoxia_heatmaps.pdf")
    }),

    tar_target(fig_numbat_heatmaps_permissive,
      # Always re-render: output layout can depend on graphics device/session state.
      make_numbat_heatmaps(
        seus_low_hypoxia, numbat_rds_files,
        p_min = 0.5, line_width = 0.1, extension = "_low_hypoxia",
        show_segment_names_on_x = TRUE, filter_midline = FALSE
      ),
      pattern = map(seus_low_hypoxia),
      iteration = "list",
      cue = tar_cue("always")  # force re-render each run: heatmap layout depends on
                                # session graphics device state, not captured by hash
    ),

    tar_target(filtered_numbat_heatmaps_file, qpdf::pdf_combine(map_chr(fig_numbat_heatmaps_permissive, 1), "results/filtered_heatmaps.pdf")),

    tar_target(filtered_large_scna_prob_file, {
      scna_var_paths <- na.omit(map_chr(fig_numbat_heatmaps_permissive, 2))
      qpdf::pdf_combine(scna_var_paths, "results/filtered_scna_probabilities.pdf")
    }),

    # --- heatmap collages ---

    tar_target(heatmap_collages,
      plot_seu_marker_heatmap_all_resolutions(seus_low_hypoxia, numbat_rds_files, large_clone_simplifications),
      pattern = map(seus_low_hypoxia),
      iteration = "list"
    ),

    tar_target(heatmap_collages_6p,
      plot_seu_marker_heatmap_all_resolutions(hypoxia_seus_6p, numbat_rds_files, large_clone_simplifications),
      pattern = map(hypoxia_seus_6p),
      iteration = "list"
    ),

    tar_target(annotated_heatmap_collages,
      plot_seu_marker_heatmap(seus_low_hypoxia, cluster_orders, numbat_rds_files, large_clone_simplifications),
      pattern = map(seus_low_hypoxia),
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
      # Persist hypoxia-scored Seurat objects for downstream thresholded analyses.
      load_and_save_hypoxia_score(filtered_seus),
      pattern = map(filtered_seus),
      iteration = "list",
      error = "null",
      cue = tar_cue(command = FALSE)
    ),

    # Per-sample threshold: reads from hypoxia_thresholds dict, defaults to 0.4.
    # Branched over hypoxia_seus so only samples with changed thresholds rebuild.
    tar_target(
      hypoxia_threshold_per_sample,
      {
        if (is.null(hypoxia_seus) || is.na(hypoxia_seus)) return(NULL)
        sample_id <- stringr::str_extract(hypoxia_seus, "SR[RX][0-9]+")
        if (!is.null(hypoxia_thresholds) && length(sample_id) == 1 && !is.na(sample_id) && sample_id %in% names(hypoxia_thresholds))
          as.numeric(hypoxia_thresholds[[sample_id]])
        else
          0.4
      },
      pattern = map(hypoxia_seus),
      iteration = "list",
      error = "null"
    ),

    tar_target(
      hypoxia_score_plots,
      plot_hypoxia_score(hypoxia_seus, threshold = hypoxia_threshold_per_sample),
      pattern = map(hypoxia_seus, hypoxia_threshold_per_sample),
      iteration = "list"
    ),

    tar_target(seus_low_hypoxia,
      # Low-hypoxia partition is later remapped into SCNA-specific subsets.
      subset_seu_by_expression(
        hypoxia_seus, run_hypoxia_clustering = TRUE,
        hypoxia_expr = glue::glue("hypoxia_score <= {hypoxia_threshold_per_sample}"),
        slug = "hypoxia_low",
        assay = "gene"
      ),
      pattern = map(hypoxia_seus, hypoxia_threshold_per_sample),
      iteration = "list",
      error = "null",
      cue = tar_cue(command = FALSE)
    ),

    tar_target(seus_high_hypoxia,
      subset_seu_by_expression(
        hypoxia_seus, run_hypoxia_clustering = TRUE,
        hypoxia_expr = glue::glue("hypoxia_score > {hypoxia_threshold_per_sample}"),
        slug = "hypoxia_high"
      ),
      pattern = map(hypoxia_seus, hypoxia_threshold_per_sample),
      iteration = "list",
      error = "null",
      cue = tar_cue(command = FALSE)
    ),

    tar_target(heatmap_collages_hypoxia,
      plot_seu_marker_heatmap(hypoxia_seus, nb_paths = numbat_rds_files, clone_simplifications = large_clone_simplifications, tmp_plot_path = TRUE),
      pattern = map(hypoxia_seus),
      iteration = "list"
    ),

    tar_target(heatmap_collages_low_hypoxia,
      plot_seu_marker_heatmap(
        seus_low_hypoxia,
        nb_paths = numbat_rds_files,
        clone_simplifications = large_clone_simplifications,
        tmp_plot_path = TRUE
      ),
      pattern = map(seus_low_hypoxia),
      iteration = "list"
    ),

    tar_target(hypoxia_gene_heatmap_low_hypoxia,
      plot_hypoxia_gene_heatmap(seus_low_hypoxia),
      pattern = map(seus_low_hypoxia),
      iteration = "list",
      error = "null"
    ),

    tar_target(heatmap_collages_high_hypoxia,
      plot_seu_marker_heatmap(
        seus_high_hypoxia,
        nb_paths = numbat_rds_files,
        clone_simplifications = large_clone_simplifications,
        tmp_plot_path = TRUE
      ),
      pattern = map(seus_high_hypoxia),
      iteration = "list"
    ),

    tar_target(hypoxia_effect_plots,
      {
        hypoxia_paths <- unlist(hypoxia_score_plots)
        hypoxia_paths <- hypoxia_paths[!is.na(hypoxia_paths)]
        heatmap_paths <- unlist(heatmap_collages_hypoxia)
        heatmap_paths <- heatmap_paths[!is.na(heatmap_paths)]
        list(
          qpdf::pdf_combine(hypoxia_paths, "01_hypoxia_score_plots.pdf"),
          qpdf::pdf_combine(heatmap_paths, "02_hypoxia_heatmap.pdf")
          # TODO: add low/high hypoxia heatmaps once those targets are fixed
        )
      },
    ),

    tar_target(silhouette_plots,
      calc_silhouette(seus_low_hypoxia),
      pattern = map(seus_low_hypoxia),
      iteration = "list"
    ),

    # --- integrated seus ---

    tar_target(
      integrated_seu_16q,
      readRDS("output/seurat/SRX11133594_SRX11133593_SRX11133592_seu.rds")
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
  # targets are named e.g. clustree_SRX10264519, clustree_SRX10264523_branch_5)
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
      # Combine static per-sample clustree branches into a single clustrees target.
      clustrees,
      .clustree_map[["clustree"]],
      command = list(!!!.x)
    )
    list(
      .clustree_map,
      .clustree_combined,
      tarchetypes::tar_files(
        clustree_compilation,
        {
          paths <- unlist(clustrees)
          paths <- paths[!is.na(paths)]
          qpdf::pdf_combine(paths, "results/clustrees.pdf")
        },
        format = "file"
      )
    )
  },

  # hypoxia_seus_1q/2p/6p/16q — subset seus_low_hypoxia to each SCNA's samples
  tarchetypes::tar_map(
    # Expand hypoxia_seus_<scna> targets from shared low-hypoxia Seurat set.
    values = scna_map_values[, "scna", drop = FALSE],
    names  = "scna",
    tar_target(hypoxia_seus,
      str_subset(unlist(seus_low_hypoxia), str_c(rb_scna_samples[[scna]], collapse = "|"))
    )
  ),

  # integrated_seu_low_hypoxia_1q/2p/6p/16q — integration within each SCNA at low hypoxia
  tarchetypes::tar_map(
    # Integrate low-hypoxia Seurat objects separately for each SCNA stratum.
    values = tibble::tibble(
      scna        = c("1q", "2p", "6p", "16q"),
      hypoxia_sym = rlang::syms(c("hypoxia_seus_1q", "hypoxia_seus_2p",
                                  "hypoxia_seus_6p", "hypoxia_seus_16q"))
    ),
    names = "scna",
    tar_target(integrated_seu_low_hypoxia,
      integration_by_scna_clones(
        hypoxia_sym, scna_of_interest = scna,
        large_clone_comparisons,
        filter_expr = "hypoxia_score <= 0.4"
      ),
      error = "null",
      cue = tar_cue(command = FALSE)
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
