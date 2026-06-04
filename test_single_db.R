library(DBI)
library(RSQLite)
library(stringr)
# library(numbatHelpers)
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")

debug(make_numbat_heatmaps)

# Query database for a sample
con <- DBI::dbConnect(RSQLite::SQLite(), 'batch_hashes.sqlite')
samples <- DBI::dbGetQuery(con,
  "SELECT DISTINCT sample_id, filepath FROM seurat_objects
   WHERE processing_stage = 'raw' AND sample_id LIKE 'SRX%' LIMIT 1")
DBI::dbDisconnect(con)

if (nrow(samples) > 0) {
  seu_path <- samples$filepath[1]
  sample_id <- samples$sample_id[1]

  cat("Testing sample:", sample_id, "\n")
  cat("Seurat path:", seu_path, "\n")

  # Find matching numbat files
  numbat_files <- c(
    list.files('output/numbat_sridhar', recursive = TRUE, pattern = '_numbat.rds$', full.names = TRUE),
    list.files('output/numbat_sridhar_filtered', recursive = TRUE, pattern = '_numbat.rds$', full.names = TRUE)
  )

  cat("Found", length(numbat_files), "numbat files\n\n")

  # Run make_numbat_heatmaps
  tryCatch({
    cat("Calling make_numbat_heatmaps...\n")
    result <- make_numbat_heatmaps(seu_path, numbat_files, extension = '_test')

    cat("\nResult:\n")
    print(result)

    for (f in result[!is.na(result)]) {
      if (file.exists(f)) {
        size_kb <- round(file.size(f) / 1024)
        cat("✓ Generated:", f, "(", size_kb, "KB)\n")
      } else {
        cat("✗ File not found:", f, "\n")
      }
    }
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    traceback()
  })
} else {
  cat("No samples found\n")
}
