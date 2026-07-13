# Audit the auto-generated cluster_dictionary for excessive cell filtering.
#
# filtered_seus drops every cell in any res.0.2 cluster the dictionary flags
# remove=1 (cluster_remove, the highest-precedence reason in
# annotate_filter_reason()). This script quantifies, per sample and per cluster,
# how many cells each remove=1 flag costs and how strong the marker evidence is,
# so we can see which samples are over-flagged.
#
# Read-only: uses the already-written discrepancies TSV (dictionary + evidence)
# and exact per-res.0.2-cluster cell counts stored in the cell_metadata table of
# batch_hashes.sqlite (column gene_snn_res.0.2 -> summary_json {cluster: n_cells}).
# No Seurat objects are loaded.

suppressMessages({
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")
  library(dplyr)
  library(readr)
  library(stringr)
})

discrep_path <- "results/cluster_dictionary/cluster_dictionary_discrepancies.tsv"
sqlite_path  <- "batch_hashes.sqlite"
out_dir      <- "doc"

stopifnot(file.exists(discrep_path))

dict <- read_tsv(discrep_path, show_col_types = FALSE) %>%
  mutate(gene_snn_res.0.2 = as.integer(`gene_snn_res.0.2`))

# ---- per-res.0.2-cluster cell counts from cell_metadata.summary_json ----------
con <- connect_hash_db(sqlite_path)
on.exit(DBI::dbDisconnect(con), add = TRUE)

counts_raw <- DBI::dbGetQuery(con,
  "SELECT filepath, summary_json
     FROM cell_metadata
    WHERE column_name = 'gene_snn_res.0.2'")

parse_counts <- function(filepath, summary_json) {
  sid <- str_extract(filepath, "SR[RX][0-9]+")
  if (is.na(sid) || is.na(summary_json) || !nzchar(summary_json)) return(NULL)
  j <- jsonlite::fromJSON(summary_json)
  tibble(
    sample_id        = sid,
    gene_snn_res.0.2 = as.integer(names(j)),
    n_cells          = as.integer(unlist(j))
  )
}

cluster_counts <- purrr::pmap_dfr(
  list(counts_raw$filepath, counts_raw$summary_json), parse_counts)

sample_totals <- cluster_counts %>%
  group_by(sample_id) %>%
  summarise(sample_total = sum(n_cells), .groups = "drop")

# ---- join counts onto the dictionary + derive evidence/confidence -------------
count_matched <- function(matched) {
  ifelse(is.na(matched) | matched == "", 0L,
         lengths(str_split(matched, ",")))
}

audit <- dict %>%
  left_join(cluster_counts, by = c("sample_id", "gene_snn_res.0.2")) %>%
  left_join(sample_totals, by = "sample_id") %>%
  mutate(
    n_matched = count_matched(matched),
    pct_of_sample = round(100 * n_cells / sample_total, 2)
  )

removed <- audit %>%
  filter(remove == "1") %>%
  mutate(
    confidence = case_when(
      abbreviation %in% c("rod", "cone") & n_matched <= 2 ~ "low",
      rule == "top1_gene"                                 ~ "low",   # incidental (e.g. APOE)
      TRUE                                                ~ "high"
    )
  ) %>%
  arrange(desc(pct_of_sample)) %>%
  transmute(
    sample_id,
    cluster = `gene_snn_res.0.2`,
    abbreviation, rule, n_matched, matched,
    n_cells, pct_of_sample, top_genes, confidence
  )

# ---- per-sample summary ------------------------------------------------------
program_levels <- c("rod", "cone", "low_qual", "MALAT1", "APOE")

by_sample <- removed %>%
  group_by(sample_id) %>%
  summarise(
    n_removed_clusters = n(),
    cells_removed = sum(n_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # bring in every sample that has counts, even with 0 removed clusters
  right_join(sample_totals, by = "sample_id") %>%
  mutate(
    n_removed_clusters = coalesce(n_removed_clusters, 0L),
    cells_removed = coalesce(cells_removed, 0L),
    pct_removed = round(100 * cells_removed / sample_total, 2)
  )

# per-program cell counts per sample
prog_wide <- removed %>%
  filter(abbreviation %in% program_levels) %>%
  group_by(sample_id, abbreviation) %>%
  summarise(cells = sum(n_cells, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = abbreviation, values_from = cells,
                     values_fill = 0, names_prefix = "cells_")

by_sample <- by_sample %>%
  left_join(prog_wide, by = "sample_id") %>%
  mutate(across(starts_with("cells_"), ~ coalesce(.x, 0L))) %>%
  arrange(desc(pct_removed))

# ---- write outputs -----------------------------------------------------------
fs::dir_create(out_dir)
write_tsv(removed,   file.path(out_dir, "cluster_dictionary_audit.tsv"), na = "NA")
write_tsv(by_sample, file.path(out_dir, "cluster_dictionary_audit_by_sample.tsv"), na = "NA")

prog_totals <- removed %>%
  count(abbreviation, name = "n_clusters") %>%
  left_join(
    removed %>% group_by(abbreviation) %>%
      summarise(cells = sum(n_cells, na.rm = TRUE), .groups = "drop"),
    by = "abbreviation"
  ) %>%
  arrange(desc(n_clusters))

low_conf <- removed %>% filter(confidence == "low")

md <- c(
  "# Cluster-dictionary filtering audit",
  "",
  sprintf("_Generated %s from `%s` + `cell_metadata` cell counts._",
          format(Sys.Date()), discrep_path),
  "",
  "Quantifies cells removed by the auto cluster_dictionary (clusters flagged",
  "`remove=1`, dropped wholesale as `cluster_remove` in `annotate_filter_reason`).",
  "",
  "## Samples ranked by % cells removed (dictionary only)",
  "",
  "| sample | total cells | clusters removed | cells removed | % removed |",
  "|---|---:|---:|---:|---:|",
  by_sample %>%
    mutate(row = sprintf("| %s | %d | %d | %d | %s |",
                         sample_id, sample_total, n_removed_clusters,
                         cells_removed, pct_removed)) %>%
    pull(row),
  "",
  "## Program totals (clusters / cells removed)",
  "",
  "| program | clusters | cells |",
  "|---|---:|---:|",
  prog_totals %>%
    mutate(row = sprintf("| %s | %d | %d |", abbreviation, n_clusters, cells)) %>%
    pull(row),
  "",
  sprintf("## Low-confidence removals (%d) — review these", nrow(low_conf)),
  "",
  "rod/cone flagged on only the 2-gene rule minimum, or incidental `top1_gene`",
  "removals (e.g. APOE). These are the likeliest false-positive deletions.",
  "",
  "| sample | cluster | abbrev | rule | n_matched | cells | % | top_genes |",
  "|---|---:|---|---|---:|---:|---:|---|",
  if (nrow(low_conf) == 0) "| _none_ | | | | | | | |" else
    low_conf %>%
      mutate(row = sprintf("| %s | %d | %s | %s | %d | %d | %s | %s |",
                           sample_id, cluster, abbreviation, rule,
                           n_matched, n_cells, pct_of_sample, top_genes)) %>%
      pull(row)
)
writeLines(md, file.path(out_dir, "cluster_dictionary_audit.md"))

# ---- console summary ---------------------------------------------------------
cat("\n=== cluster_dictionary audit ===\n")
cat(sprintf("samples: %d | total remove=1 clusters: %d | low-confidence: %d\n",
            nrow(by_sample), nrow(removed), nrow(low_conf)))
cat("\nTop over-flagged samples (by %% removed):\n")
print(head(by_sample, 12), n = 12)
cat("\nProgram totals:\n"); print(prog_totals)
cat("\nLow-confidence removals:\n")
print(low_conf, n = nrow(low_conf))
cat(sprintf("\nWrote:\n  %s/cluster_dictionary_audit.tsv\n  %s/cluster_dictionary_audit_by_sample.tsv\n  %s/cluster_dictionary_audit.md\n",
            out_dir, out_dir, out_dir))
