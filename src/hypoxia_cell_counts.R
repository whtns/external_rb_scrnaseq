#!/usr/bin/env Rscript
# Compile per-sample low/high hypoxia cell counts from the barcodes DB.
# Does not load any Seurat RDS files.
# Output: results/hypoxia_cluster_split/hypoxia_cell_counts.csv

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(tidyr)
  library(readr)
})

con <- dbConnect(SQLite(), "batch_hashes.sqlite")

rows <- dbGetQuery(con,
  "SELECT sample_id, seu_type, cells
   FROM seu_cells
   WHERE seu_type IN ('hypoxia_low', 'hypoxia_high')
   ORDER BY sample_id, seu_type")

dbDisconnect(con)

rows$n_cells <- vapply(rows$cells, function(x) {
  length(strsplit(x, "\n", fixed = TRUE)[[1]])
}, integer(1L))

counts <- rows |>
  select(sample_id, seu_type, n_cells) |>
  pivot_wider(names_from = seu_type, values_from = n_cells,
              values_fill = NA_integer_) |>
  rename(n_low = hypoxia_low, n_high = hypoxia_high) |>
  mutate(
    n_total  = rowSums(cbind(n_low, n_high), na.rm = TRUE),
    pct_high = round(100 * n_high / n_total, 1)
  ) |>
  arrange(sample_id)

print(counts, n = Inf)

out <- "results/hypoxia_cluster_split/hypoxia_cell_counts.csv"
write_csv(counts, out)
cat("\nWritten to", out, "\n")
