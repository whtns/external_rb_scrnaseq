library(targets)
library(stringr)
library(dplyr)
library(Seurat)
library(numbatHelpers)

# Load one sample's data from targets
tar_load(unfiltered_seus, store = "_targets_r431")
tar_load(numbat_rds_files, store = "_targets_r431")

cat("Loaded unfiltered_seus, length:", length(unfiltered_seus), "\n")
cat("Loaded numbat_rds_files, length:", length(numbat_rds_files), "\n")

# Test make_numbat_heatmaps on the first sample
if (length(unfiltered_seus) > 0 && length(numbat_rds_files) > 0) {
  seu_path <- unfiltered_seus[[1]]
  cat("\nTesting make_numbat_heatmaps with seu_path:", seu_path, "\n")
  
  tryCatch({
    result <- make_numbat_heatmaps(
      seu_path,
      numbat_rds_files,
      p_min = 0.9,
      line_width = 0.1,
      extension = "_unfiltered_test",
      show_segment_names_on_x = TRUE
    )
    cat("Success! Result:\n")
    print(result)
    cat("\nResult class:", class(result), "\n")
    cat("Result length:", length(result), "\n")
    cat("Result[1]:", result[1], "\n")
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    traceback()
  })
} else {
  cat("Missing data\n")
}
