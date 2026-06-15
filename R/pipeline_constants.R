# Pipeline-level constants used across pipeline_targets_*.R files.
# Defined outside tar_plan() so they are available at pipeline-definition time.

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

# Values tibble for tarchetypes::tar_map() over the 4 SCNA types.
# Columns:
#   scna      - SCNA label used in target names and as list index key
#   seus_sym  - symbol of the corresponding low-hypoxia Seurat path target
#   var_y     - y-axis grouping for clone_pearls plots
scna_map_values <- tibble::tibble(
  scna     = c("1q",         "2p",       "6p",       "16q"),
  seus_sym = rlang::syms(c(
    "hypoxia_seus_1q", "hypoxia_seus_2p",
    "hypoxia_seus_6p", "hypoxia_seus_16q"
  )),
  var_y    = c("phase_level", "clusters", "clusters", "phase_level")
)

# rclone destination for figures_and_tables sync.
# Configure with: module load rclone && rclone config
gdrive_destination <- "gdrive:rb_scrnaseq/figures_and_tables"

# Document-order mapping: semantic target name → display label.
# Edit this vector to renumber figures without touching target definitions.
figure_order <- c(
  # Main figures
  fig_single_sample_panels               = "Fig. 2",
  fig_16q_alternative_resolutions        = "Fig. 3",
  fig_1q_integrated                      = "Fig. 4",
  fig_6p_integrated                      = "Fig. 5",
  fig_1q_cluster_diffex_integrated       = "Fig. 7a",
  fig_1q_cluster_diffex_unintegrated     = "Fig. 7b",
  fig_16q_cluster_diffex_integrated      = "Fig. 8a",
  fig_16q_cluster_diffex_unintegrated    = "Fig. 8b",
  fig_2p_corresponding_clusters          = "Fig. 9",
  fig_6p_corresponding_clusters          = "Fig. 10",
  # Supplemental figures
  fig_tcga_scna_frequency                = "Fig. S2",
  fig_tcga_gistic                        = "Fig. S3",
  fig_numbat_heatmaps                    = "Fig. S3a",
  fig_study_cell_stats                   = "Fig. S4",
  fig_numbat_expression_smoothed         = "Fig. S5",
  fig_karyograms                         = "Fig. S6a",
  fig_1q_sample_specific_integrated      = "Fig. S7",
  fig_1q_clone_diffex_within_clusters    = "Fig. S8",
  fig_1q_cluster_diffex_of_interest      = "Fig. S9",
  fig_16q_clone_diffex_within_clusters   = "Fig. S10",
  fig_16q_sample_specific_integrated     = "Fig. S12",
  fig_numbat_heatmaps_permissive         = "Fig. S13",
  fig_2p_clone_diffex_within_clusters    = "Fig. S20",
  fig_6p_sample_specific_integrated      = "Fig. S23",
  fig_subtype_markers                    = "Fig. S25",
  fig_2p_sample_specific_integrated      = "Fig. S4.8",
  fig_2p_sample_specific_unintegrated    = "Fig. S4.9",
  fig_1q_integrated_v2                   = "Fig. 4v2",
  fig_1q_16q_combined                    = "Fig. 4-7",
  fig_2p_integrated                      = "Fig. 4-9",
  fig_regression_diagnostics             = "Fig. 3-5",
  fig_single_sample_panels_with_diploid  = "Fig. 2 (diploid)",
  # Tables
  table_qc_stats                         = "Table S1",
  table_sample_metadata                  = "Table S2",
  table_removed_clusters                 = "Table S3",
  table_rod_rich_samples                 = "Table 6",
  table_1q_clone_per_cluster             = "Table S7",
  table_16q_clone_per_cluster            = "Table S9",
  table_2p_clone_per_cluster             = "Table S10"
)

# Values tibble for tarchetypes::tar_map() over debranched sample/branch IDs.
# Columns:
#   id       - sample/branch label used in target names (e.g. clustree_SRX10264519)
#   seu_path - path to the filtered Seurat RDS for that sample/branch
debranched_map_values <- tibble::tibble(
  id = c(
    "SRX10264519",          "SRX10264520",
    "SRX10264523_branch_5", "SRX10264523_branch_6",
    "SRX10264524_branch_6", "SRX10264524_branch_4", "SRX10264524_branch_5",
    "SRX10264525",          "SRX10264526",
    "SRX11133594",          "SRX11133593",          "SRX11133592",
    "SRX11133588_branch_2", "SRX11133588_branch_3",
    "SRX11133587_branch_4", "SRX11133587_branch_7",
    "SRX11133585_branch_3", "SRX11133585_branch_4",
    "SRX14116947",          "SRX14116944",
    "SRX22868105",
    "SRX22868102_branch_3", "SRX22868102_branch_4"
  )
) |>
  dplyr::mutate(seu_path = paste0("output/seurat/", id, "_filtered_seu.rds"))
