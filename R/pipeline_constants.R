# Pipeline-level constants used across pipeline_targets_*.R files.
# Defined outside tar_plan() so they are available at pipeline-definition time.

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

hypoxia_threshold <- 0.5

# Values tibble for tarchetypes::tar_map() over the 4 SCNA types.
# Columns:
#   scna      - SCNA label used in target names and as list index key
#   seus_sym  - symbol of the corresponding debranched Seurat path target
#   var_y     - y-axis grouping for clone_pearls plots
scna_map_values <- tibble::tibble(
  scna     = c("1q",         "2p",       "6p",       "16q"),
  seus_sym = rlang::syms(c(
    "debranched_seus_1q", "debranched_seus_2p",
    "debranched_seus_6p", "debranched_seus_16q"
  )),
  var_y    = c("phase_level", "clusters", "clusters", "phase_level")
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
