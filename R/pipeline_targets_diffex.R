# Differential expression targets: clone diffex (all/cis/trans), clustree diffex,
# and integration cluster diffex (fig 07/08 inputs and plots).
# Defines: pipeline_targets_diffex (spliced into tar_plan in _targets.R)

pipeline_targets_diffex <- list(

  # --- all locations ---

  tar_target(all_diffex_clones,
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "all"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  tar_target(all_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "all"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  tar_target(
    volcano_thresholds_all,
    yaml::read_yaml(here::here("config/volcano_thresholds_all.yaml"))
  ),

  tar_target(
    volcano_all,
    make_volcano_diffex_clones(
      all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_all.pdf",
      all_diffex_clones,
      "results/diffex_bw_clones_all.pdf"
    ),
  ),

  tar_target(
    table_all_diffex_clones,
    tabulate_diffex_clones(
      all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_all.xlsx",
      "results/diffex_bw_clones_per_cluster_all_by_chr.xlsx",
      all_diffex_clones,
      "results/diffex_bw_clones_all.xlsx",
      "results/diffex_bw_clones_all_by_chr.xlsx"
    ),
  ),

  # --- cis (in-segment) ---

  tar_target(cis_diffex_clones,
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "cis"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  tar_target(cis_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "cis"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  tar_target(
    volcano_thresholds_in_segment,
    yaml::read_yaml(here::here("config/volcano_thresholds_in_segment.yaml"))
  ),

  tar_target(
    volcano_cis,
    make_volcano_diffex_clones(
      cis_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
      cis_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.pdf"
    ),
  ),

  tar_target(
    table_cis_diffex_clones,
    tabulate_diffex_clones(
      cis_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.xlsx",
      "results/diffex_bw_clones_per_cluster_large_in_segment_by_chr.xlsx",
      cis_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.xlsx",
      "results/diffex_bw_clones_large_in_segment_by_chr.xlsx"
    ),
  ),

  # --- trans (out-of-segment) ---

  tar_target(trans_diffex_clones,
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "out_of_segment"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  tar_target(trans_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "out_of_segment"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  tar_target(
    volcano_trans,
    make_volcano_diffex_clones(
      trans_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf",
      trans_diffex_clones,
      "results/diffex_bw_clones_large_out_of_segment.pdf"
    ),
  ),

  tar_target(
    found_kooi_candidates,
    tally_kooi_candidates(
      cis_diffex_clones  = "results/diffex_bw_clones_large_in_segment_by_chr.xlsx",
      trans_diffex_clones = "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx"
    )
  ),

  tar_target(
    table_trans_diffex_clones,
    tabulate_diffex_clones(
      trans_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_per_cluster_large_out_of_segment_by_chr.xlsx",
      trans_diffex_clones,
      "results/diffex_bw_clones_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx"
    ),
  ),

  # --- clustree-based diffex ---

  tar_target(
    clustree_tables,
    pull_clustree_tables(clustrees, divergent_cluster_file)
  ),

  tar_target(clustree_diffexes,
    find_all_diffex_from_clustree(clustree_tables, debranched_seus, clone_comparisons = large_clone_comparisons),
    pattern = map(clustree_tables),
  ),

  tar_target(
    clustree_cis_changes,
    find_candidate_cis_in_clustree_diffexes(clustree_diffexes)
  ),

  tar_target(
    clustree_trans_changes,
    find_candidate_trans_in_clustree_diffexes(clustree_diffexes)
  ),

  tar_target(
    clustree_all_changes,
    find_candidate_all_in_clustree_diffexes(clustree_diffexes)
  ),

  tar_target(
    num_diffex_clone_scna_tally,
    tally_num_diffex(unfiltered_oncoprint_input_by_scna)
  ),

  # --- integration cluster diffex: 1q (fig 07) ---

  tar_target(fig_07a_input,  # 1q+ cluster diffex after integration
    find_diffex_bw_clones_for_each_cluster(integrated_seus_1q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "1q+"),
    pattern = map(integrated_seus_1q),
    iteration = "list"
  ),

  tar_target(fig_07a,
    plot_fig_07_08(fig_07a_input, plot_path = "results/fig_07a.pdf", plot_title = "fig_07a: 1q+ cluster diffex after integration", height = 5, width = 4),
    iteration = "list"
  ),

  tar_target(fig_07b_input,  # 1q+ cluster diffex without integration
    find_diffex_bw_clones_for_each_cluster(debranched_seus_1q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "1q+"),
    pattern = map(debranched_seus_1q),
    iteration = "list"
  ),

  tar_target(fig_07b,
    plot_fig_07_08(fig_07b_input, plot_path = "results/fig_07b.pdf", plot_title = "fig_07b: 1q+ cluster diffex without integration", height = 5, width = 4),
    iteration = "list"
  ),

  # --- integration cluster diffex: 16q (fig 08) ---

  tar_target(fig_08a_input,  # 16q- cluster diffex after integration
    find_diffex_bw_clones_for_each_cluster(integrated_seus_16q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-"),
    pattern = map(integrated_seus_16q),
    iteration = "list"
  ),

  tar_target(fig_08a,
    plot_fig_07_08(fig_08a_input, plot_title = "fig_08a: 16q- cluster diffex after integration", plot_path = "results/fig_08a.pdf", height = 5, width = 6),
    iteration = "list"
  ),

  tar_target(fig_08b_input,  # 16q- cluster diffex without integration
    find_diffex_bw_clones_for_each_cluster(debranched_seus_16q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-"),
    pattern = map(debranched_seus_16q),
    iteration = "list"
  ),

  tar_target(fig_08b,
    plot_fig_07_08(fig_08b_input, plot_title = "fig_08b: 16q- cluster diffex without integration", plot_path = "results/fig_08b.pdf", p_adj_threshold = 0.05, height = 5, width = 6),
    iteration = "list"
  )

)
