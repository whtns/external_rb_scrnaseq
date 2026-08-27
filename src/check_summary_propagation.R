# Dry run: does the majority-round numbat RDS change actually reach sample_summaries?
#
# Two questions:
#   1. Does numbat_rds_files see the rebuild? (it is tar_files/format="file" and
#      freeze_rds_files is FALSE, so it should)
#   2. Do the depend = FALSE cues on unfiltered_seus / filtered_seus stop that
#      change from propagating downstream to sample_summaries?
#
# Read-only. Invalidates nothing, builds nothing.
suppressPackageStartupMessages({
  library(targets)
  library(data.table)
})
tar_config_set(store = "_targets_r431")

cat("=== 1. RDS mtime vs what the store last recorded ===\n")
m  <- as.data.table(tar_meta())
nb <- m[grepl("^numbat_rds_files", name) & !is.na(time)]
cat("numbat_rds_files branches recorded:", nrow(nb), "\n")
if (nrow(nb)) {
  paths <- vapply(nb$path, function(p) if (length(p)) p[1] else NA_character_, "")
  info  <- file.info(paths)
  newer <- which(!is.na(info$mtime) & info$mtime > nb$time)
  cat("RDS modified AFTER the store recorded them:", length(newer), "\n")
  if (length(newer)) print(basename(paths[newer]))
}

cat("\n=== 2. outdated upstream of sample_summaries (authoritative) ===\n")
od <- tar_outdated(names = any_of("sample_summaries"), reporter = "silent")
cat("outdated targets:", length(od), "\n")
print(sort(od))

cat("\n=== 3. are the frozen seu targets in that outdated set? ===\n")
for (t in c("numbat_rds_files", "unfiltered_seus", "filtered_seus",
            "filtered_seus_with_phase", "hypoxia_partition_paths",
            "seus_low_hypoxia", "unfiltered_clone_tree_files",
            "numbat_heatmap_plots_unfiltered", "sample_summaries")) {
  cat(sprintf("  %-32s outdated: %s\n", t, t %in% od))
}
cat("\nIf numbat_rds_files is outdated but unfiltered_seus is NOT,\n",
    "the depend = FALSE cue is the reason and it needs an explicit tar_invalidate.\n", sep = "")
