library(numbatHelpers)
library(stringr)

sqlite_path <- "batch_hashes.sqlite"

filtered_paths <- fs::dir_ls("output/seurat", regexp = "_filtered_seu\\.rds$")
hypoxia_paths  <- fs::dir_ls("output/seurat", regexp = "_(hypoxia_low|hypoxia_high)_seu\\.rds$")

backfill <- function(paths, seu_type_fn) {
  for (p in paths) {
    sample_id <- str_extract(p, "SR[RX][0-9]+")
    seu_type  <- seu_type_fn(p)
    seu <- tryCatch(readRDS(p), error = function(e) {
      cat("SKIP (unreadable):", as.character(p), "-", conditionMessage(e), "\n")
      NULL
    })
    if (is.null(seu)) next
    cells <- colnames(seu)
    cat("Saving", sample_id, "/", seu_type, "->", length(cells), "cells\n")
    save_cell_barcodes_to_db(as.character(p), sample_id, seu_type, cells, sqlite_path)
  }
}

backfill(filtered_paths, function(p) "filtered")
backfill(hypoxia_paths,  function(p) str_extract(p, "hypoxia_low|hypoxia_high"))

# Verify
con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
res <- DBI::dbGetQuery(con, "SELECT sample_id, seu_type, length(cells) - length(replace(cells, char(10), '')) + 1 AS n_cells FROM seu_cells ORDER BY sample_id, seu_type")
DBI::dbDisconnect(con)
print(res)
