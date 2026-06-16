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
  # Main figures (per figure_and_table_captions.txt)
  fig_single_sample_panels               = "Fig. 2",
  # Fig. 3 = cluster marker gene analysis in single pilot tumor — covered by fig_01 (partial)
  fig_1q_integrated                      = "Fig. 4",
  fig_16q_integrated                     = "Fig. 5",
  fig_1q_cluster_diffex_integrated       = "Fig. 7a",
  fig_1q_cluster_diffex_unintegrated     = "Fig. 7b",
  fig_16q_cluster_diffex_integrated      = "Fig. 8a",
  fig_16q_cluster_diffex_unintegrated    = "Fig. 8b",
  fig_2p_corresponding_clusters          = "Fig. 9",
  fig_6p_corresponding_clusters          = "Fig. 10",
  # Supplemental figures (document order per figure_and_table_captions.txt)
  fig_tcga_scna_frequency                = "Fig. S1",
  fig_tcga_gistic                        = "Fig. S2",
  fig_numbat_heatmaps                    = "Fig. S3",
  fig_numbat_expression_smoothed         = "Fig. S3b",      # sub-panel of S3; no direct manuscript cite
  fig_study_cell_stats                   = "Fig. S4",
  fig_1q_sample_specific_integrated      = "Fig. S5",
  fig_karyograms                         = "Fig. S6a",
  # Fig. S6  = 1q+ sample-specific without integration — no pipeline target
  # Fig. S7  = alt Louvain resolutions for integrated 1q+ — no pipeline target
  fig_1q_clone_diffex_within_clusters    = "Fig. S8",
  fig_1q_cluster_diffex_of_interest      = "Fig. S9",
  fig_16q_clone_diffex_within_clusters   = "Fig. S10",
  # Fig. S11 = alt resolutions for integrated 16q- — no pipeline target
  fig_16q_sample_specific_integrated     = "Fig. S12",
  # Fig. S13 = 16q- sample-specific without integration — no pipeline target
  fig_2p_integrated                      = "Fig. S14",
  # Fig. S15 = alt resolutions for integrated 2p+ — no pipeline target
  fig_2p_sample_specific_integrated      = "Fig. S16",
  # Fig. S17 = 2p+ cluster DE — no pipeline target
  # Fig. S18 = 2p+ candidate drivers — no pipeline target
  # Fig. S19 = 6p+ SCNA boundaries — no pipeline target
  # Fig. S20 = 6p+ sample-specific without integration — no pipeline target
  # Fig. S21 = 6p+ enriched terms — no pipeline target
  # Fig. S22 = 6p+ DE in cis — no pipeline target
  fig_subtype_markers                    = "Fig. S25",
  # Draft / unassigned pipeline figures (document position not yet confirmed)
  fig_6p_integrated                      = "Fig. 6p-draft",   # was Fig. 5; captions Fig. 5 = 16q-
  fig_numbat_heatmaps_permissive         = "Fig. S-nbt-perm", # was Fig. S13; S13 = 16q- unintegrated
  fig_2p_clone_diffex_within_clusters    = "Fig. S-2p-clde",  # was Fig. S20; S20 = 6p+ unintegrated
  fig_6p_sample_specific_integrated      = "Fig. S-6p-si",    # no confirmed position in captions
  clone_tree_collage                     = "Fig. 1b",   # static PNG; also cited as Fig. S24
  clone_tree_collage_of_merged_replicates = "Fig. 1b-rep", # static PNG (not yet created)
  fig_2p_sample_specific_unintegrated    = "Fig. S4.9",
  fig_1q_integrated_v2                   = "Fig. 4v2",
  fig_1q_16q_combined                    = "Fig. 4-7",
  fig_regression_diagnostics             = "Fig. 3-5",
  fig_single_sample_panels_with_diploid  = "Fig. 2 (diploid)",
  # Tables (per figure_and_table_captions.txt)
  table_sample_metadata                  = "Table S2",
  table_rod_rich_samples                 = "Table S3",          # S3 = rod cell proportions
  table_qc_stats                         = "Table S4",          # S4 = Tumor QC
  table_removed_clusters                 = "Table S8",          # S8 = tally of removed clusters
  table_1q_clone_per_cluster             = "Table X2",          # X2 = clone per cluster (1q+ subtable A)
  table_16q_clone_per_cluster            = "Table X2 (16q-)",   # X2 subtable B
  table_2p_clone_per_cluster             = "Table X2 (2p+)"     # X2 subtable C
)

# Values tibble for tarchetypes::tar_map() over debranched sample/branch IDs.
# Columns:
#   id       - sample/branch label used in target names (e.g. clustree_SRX10264519)
#   seu_path - path to the filtered Seurat RDS for that sample/branch
debranched_map_values <- tibble::tibble(
  id = c(
    "SRX10264519",          "SRX10264520",
    "SRX10264523",
    "SRX10264524",
    "SRX10264525",          "SRX10264526",
    "SRX11133594",          "SRX11133593",          "SRX11133592",
    "SRX11133588",
    "SRX11133587",
    "SRX11133585",
    "SRX14116947",          "SRX14116944",
    "SRX22868105",
    "SRX22868102"
  )
) |>
  dplyr::mutate(seu_path = paste0("output/seurat/", id, "_filtered_seu.rds"))
