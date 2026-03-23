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
#   id       - sample/branch label used in target names (e.g. clustree_SRR13884242)
#   seu_path - path to the filtered Seurat RDS for that sample/branch
debranched_map_values <- tibble::tibble(
  id = c(
    "SRR13884242",          "SRR13884243",
    "SRR13884246_branch_5", "SRR13884246_branch_6",
    "SRR13884247_branch_6", "SRR13884247_branch_4", "SRR13884247_branch_5",
    "SRR13884248",          "SRR13884249",
    "SRR14800534",          "SRR14800535",          "SRR14800536",
    "SRR14800540_branch_2", "SRR14800540_branch_3",
    "SRR14800541_branch_4", "SRR14800541_branch_7",
    "SRR14800543_branch_3", "SRR14800543_branch_4",
    "SRR17960481",          "SRR17960484",
    "SRR27187899",
    "SRR27187902_branch_3", "SRR27187902_branch_4"
  )
) |>
  dplyr::mutate(seu_path = paste0("output/seurat/", id, "_filtered_seu.rds"))
