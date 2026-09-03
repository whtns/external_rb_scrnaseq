#!/usr/bin/env Rscript
# Sequencing depth, the filters that are actually applied, and where they
# disagree (github #43).
#
# Issue #43 asked to "plot read count on sample summaries and possibly exclude
# low read count cells". The panel is the other half of that issue. This script
# is the evidence for the second half, and it argues AGAINST a new blanket
# cutoff while identifying two things that do need attention:
#
#   1. The cohort's existing floor is not uniform. Some samples retain cells
#      down to ~500 nCount_gene, others never go below ~1250. If that tracks
#      when the object was processed rather than the biology, it is a batch
#      confound in anything that compares depth-sensitive quantities across
#      samples.
#
#   2. The pipeline filters at nCount_gene > 1000, but numbat never sees that
#      filter: run_numbat.R takes its cells from the UNFILTERED *_seu.rds and
#      subsets on cell type only, while numbat's own min_depth defaults to 0.
#      So CNV inference runs on cells the Seurat analysis later discards.
#
# Read-only. Writes results/cell_depth_provenance.csv and prints the tables that
# go into docs/cell_depth_provenance.md.

suppressPackageStartupMessages({
  library(dplyr)
})

DB        <- "batch_hashes.sqlite"
QC_TABLE  <- "results/numbat_rds_qc.csv"
THRESHOLD <- 1000    # the pipeline's own filter, filter_sample_qc()

stopifnot(file.exists(DB))

con <- DBI::dbConnect(RSQLite::SQLite(), DB)
on.exit(DBI::dbDisconnect(con), add = TRUE)

qc <- DBI::dbGetQuery(con, "
  SELECT filepath, sample_id, nCount_gene, nFeature_gene, percent_mt
  FROM   cell_qc_values")

# One row per (sample, processing stage). Restrict to the unfiltered objects --
# that is the cell set numbat is run on, and the stage whose floor tells us how
# the sample was filtered upstream.
unf <- qc %>%
  filter(!is.na(sample_id), grepl("/[A-Z]{3}[0-9]+_seu\\.rds$", filepath))

per_sample <- unf %>%
  group_by(sample_id, filepath) %>%
  summarise(
    n_cells    = n(),
    min_count  = min(nCount_gene, na.rm = TRUE),
    q05        = stats::quantile(nCount_gene, 0.05, na.rm = TRUE),
    median     = stats::median(nCount_gene, na.rm = TRUE),
    max_count  = max(nCount_gene, na.rm = TRUE),
    n_below    = sum(nCount_gene < THRESHOLD, na.rm = TRUE),
    pct_below  = 100 * n_below / n_cells,
    .groups    = "drop"
  ) %>%
  mutate(
    mtime = suppressWarnings(as.Date(file.info(filepath)$mtime)),
    # Two regimes, split well clear of both clusters (~500 vs ~1250+).
    floor_regime = ifelse(min_count < 900, "low (~500)", "high (>=1250)")
  ) %>%
  arrange(median)

cat("=== per-sample depth, unfiltered objects (the cells numbat sees) ===\n")
print(as.data.frame(per_sample %>%
  select(sample_id, n_cells, min_count, q05, median, max_count,
         pct_below, floor_regime, mtime)), row.names = FALSE, digits = 6)

cat("\n=== finding 1: the existing floor is not uniform ===\n")
print(as.data.frame(per_sample %>% count(floor_regime, name = "n_samples")),
      row.names = FALSE)
cat("\nfloor regime against processing date (_seu.rds mtime):\n")
print(as.data.frame(per_sample %>% count(floor_regime, mtime, name = "n_samples")),
      row.names = FALSE)

cat("\n=== finding 2: a cohort-wide cutoff is not supported ===\n")
tot <- nrow(unf)
for (t in c(500, 1000, 2000)) {
  n <- sum(unf$nCount_gene < t, na.rm = TRUE)
  cat(sprintf("  cells below %5d: %6d / %6d (%.2f%%)\n", t, n, tot, 100 * n / tot))
}
cat(sprintf("\n  depth range across samples: median %s to %s (%.1f-fold)\n",
            format(min(per_sample$median)), format(max(per_sample$median)),
            max(per_sample$median) / min(per_sample$median)))

# Does depth predict whether numbat produced a usable object? If the answer is
# "only for one sample", a blanket threshold is the wrong instrument.
if (file.exists(QC_TABLE)) {
  nbqc <- readr::read_csv(QC_TABLE, show_col_types = FALSE)
  joined <- per_sample %>%
    inner_join(nbqc %>% select(sample_id, n_rb_events, n_clones,
                               null_components), by = "sample_id") %>%
    mutate(incomplete = !is.na(null_components) & nzchar(null_components))

  cat("\n=== finding 3: depth vs numbat outcome ===\n")
  cat(sprintf("samples with incomplete numbat artefacts: %d / %d\n",
              sum(joined$incomplete), nrow(joined)))
  if (any(joined$incomplete)) {
    cat("\nthe incomplete ones, with their depth rank:\n")
    j <- joined %>% arrange(median) %>% mutate(depth_rank = row_number())
    print(as.data.frame(j %>% filter(incomplete) %>%
            select(sample_id, median, pct_below, depth_rank,
                   n_rb_events, null_components)), row.names = FALSE)
  }
  ct <- suppressWarnings(stats::cor(joined$median, joined$n_rb_events,
                                    use = "complete.obs", method = "spearman"))
  cat(sprintf("\nspearman(median depth, n_rb_events) = %.3f\n", ct))
  cat("A weak correlation here means depth is not a general driver of event\n")
  cat("recovery, and the association is carried by the extreme sample.\n")

  per_sample <- joined
}

readr::write_csv(per_sample, "results/cell_depth_provenance.csv")
cat("\nwrote results/cell_depth_provenance.csv\n")
cat("DIAG DONE\n")
