library(stringr)
library(numbatHelpers)

# Find actual Seurat files
seu_files <- list.files("output/seurat", pattern = ".rds$", full.names = TRUE)
cat("Found", length(seu_files), "Seurat files\n")

# Find numbat files
numbat_files <- c(
  list.files("output/numbat", recursive = TRUE, pattern = "_numbat.rds$", full.names = TRUE),
  list.files("output/numbat_sridhar", recursive = TRUE, pattern = "_numbat.rds$", full.names = TRUE),
  list.files("output/numbat_sridhar_filtered", recursive = TRUE, pattern = "_numbat.rds$", full.names = TRUE)
)
cat("Found", length(numbat_files), "numbat files\n")

if (length(seu_files) > 0 && length(numbat_files) > 0) {
  # Use first Seurat file
  seu_path <- seu_files[1]
  sample_id <- str_extract(seu_path, "SRX[0-9]+")
  cat("\nTesting with seu_path:", seu_path, "\n")
  cat("Sample ID:", sample_id, "\n")

  # Create test directory
  dir.create("results/clone_tree_test", showWarnings = FALSE)

  # Run make_numbat_heatmaps
  tryCatch({
    cat("\nCalling make_numbat_heatmaps...\n")
    result <- make_numbat_heatmaps(
      seu_path,
      numbat_files,
      p_min = 0.9,
      line_width = 0.1,
      extension = "_test",
      show_segment_names_on_x = TRUE
    )
    cat("SUCCESS! Result:\n")
    print(result)

    # Check if files were created
    for (file in result) {
      if (!is.na(file)) {
        cat("Checking file:", file, "\n")
        if (file.exists(file)) {
          cat("  ✓ File exists (", file.size(file), "bytes)\n")
        } else {
          cat("  ✗ File NOT found\n")
        }
      }
    }
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    traceback()
  })
} else {
  cat("Missing test data\n")
}
