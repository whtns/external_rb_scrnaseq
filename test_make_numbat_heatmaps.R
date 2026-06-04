library(stringr)
library(dplyr)
library(Seurat)
library(numbatHelpers)

# Find actual data files
unfiltered_seus_files <- Sys.glob("_targets_r431/objects/unfiltered_seus_*")
unfiltered_seus_file <- head(unfiltered_seus_files[!grepl("/", unfiltered_seus_files)], 1)

numbat_rds_files <- Sys.glob("data/numbat_output/*/numbat_object.rds")

cat("Testing make_numbat_heatmaps with:\n")
cat("  unfiltered seu:", unfiltered_seus_file, "\n")
cat("  numbat rds files:", length(numbat_rds_files), "files\n")

if (length(numbat_rds_files) == 0) {
  cat("No numbat RDS files found. Looking in alternative locations...\n")
  numbat_rds_files <- Sys.glob("data/numbat_*.rds")
  cat("  Found:", length(numbat_rds_files), "files\n")
}

if (length(unfiltered_seus_file) == 0 || length(numbat_rds_files) == 0) {
  cat("Missing test data files.\n")
  quit("no")
}

unfiltered_seu_path <- unfiltered_seus_file[1]

tryCatch({
  result <- make_numbat_heatmaps(
    unfiltered_seu_path,
    numbat_rds_files,
    p_min = 0.9,
    line_width = 0.1,
    extension = "_unfiltered",
    show_segment_names_on_x = TRUE
  )
  cat("Success! Result:\n")
  print(result)
  cat("Result class:", class(result), "\n")
  cat("Result length:", length(result), "\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})
