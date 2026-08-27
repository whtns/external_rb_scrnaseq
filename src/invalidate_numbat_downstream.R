# Free every target that consumes numbat RDS *content* but is gated from seeing
# a content change. Run before rebuilding sample_summaries after a numbat RDS
# migration.
#
# Two independent gates make the store lie about freshness:
#
#   1. cue = tar_cue(depend = FALSE) on unfiltered_seus / filtered_seus.
#      tar_outdated() reports these as outdated, but tar_make() applies the cue
#      at build time and dispatches zero branches. Observed 2026-08-26: both
#      patterns were declared, no branch ran, metadata stayed at 2026-06-30.
#
#   2. Path-string inputs. nb_paths_s06a is dir_ls() returning a character
#      vector, NOT format = "file", so identical filenames hash identically no
#      matter how the RDS content changed -- ideogram_res_s06a_multi never
#      re-runs. The same applies to every target whose upstream passes a saved
#      path (hypoxia_partition_paths, seus_low_hypoxia, the clone-tree and
#      heatmap targets, sample_summaries itself).
#
# Invalidation drops metadata only; no data file is deleted. Targets rebuild and
# overwrite their outputs on the next tar_make().
suppressPackageStartupMessages(library(targets))
tar_config_set(store = "_targets_r431")

TARGETS <- c(
  # gate 1: depend = FALSE
  "unfiltered_seus", "filtered_seus", "filter_inspection_metadata",
  "filtered_seus_with_phase",

  # hypoxia chain (passes saved paths at every hop)
  "hypoxia_seus", "hypoxia_threshold_per_sample", "hypoxia_partition_paths",
  "seus_low_hypoxia", "seus_high_hypoxia",

  # gate 2: ideograms, path-gated behind nb_paths_s06a
  "nb_paths_s06a", "ideogram_res_s06a_multi", "ideogram_res_s06a_unfiltered",
  "ideogram_res_s06a_filtered", "ideogram_res_s06a_low_hypoxia",

  # clone trees / heatmaps / expression: rebuilt off numbat_rds_files already,
  # but they also take the seu objects, which are about to change
  "unfiltered_clone_tree_files", "unfiltered_clone_trees_segments_files",
  "filtered_clone_tree_files", "filtered_clone_trees_segments_files",
  "low_hypoxia_clone_tree_files", "low_hypoxia_clone_trees_segments_files",
  "numbat_heatmap_plots_unfiltered", "numbat_heatmap_plots_subset",
  "numbat_heatmap_plots_low_hypoxia",
  "low_hypoxia_numbat_expression", "low_hypoxia_numbat_bulk_clones",
  "filtering_cell_counts_table",

  # the deliverable
  "sample_summaries"
)

meta <- tar_meta(complete_only = FALSE)
meta$parent <- sub("_[0-9a-f]{16}$", "", meta$name)

cat("=== before ===\n")
for (t in TARGETS) {
  s <- meta[meta$parent == t & !is.na(meta$time), ]
  cat(sprintf("  %-40s branches=%3d latest=%s\n", t, nrow(s),
              if (nrow(s)) as.character(max(s$time)) else "absent"))
}

tar_invalidate(names = any_of(TARGETS))

after <- tar_meta(complete_only = FALSE)
after$parent <- sub("_[0-9a-f]{16}$", "", after$name)
left <- sum(after$parent %in% TARGETS)
cat("\ninvalidated", length(TARGETS), "target names;",
    "metadata records remaining for them:", left, "\n")
