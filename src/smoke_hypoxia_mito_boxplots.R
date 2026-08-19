#!/usr/bin/env Rscript
# Smoke test for plot_hypoxia_mito_score_boxplots() on three samples.
#
# The point of the test is the RECONCILIATION, not just "a PDF appeared": the
# replay re-clusters the hypoxia objects itself, so the only evidence that its
# clusters are the split's clusters is that n_cells and mean_score match the
# logged values exactly. A 100% reconciliation rate on the joined rows is the
# pass condition; anything less means the replay drifted.

suppressPackageStartupMessages({
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")
  library(dplyr)
})

log_csv <- "results/hypoxia_cluster_split/hypoxia_split_log_all.csv"
stopifnot(file.exists(log_csv))

seus <- sort(Sys.glob("output/seurat/*_seu_hypoxia.rds"))
seus <- seus[!grepl("_labeled", seus)]
cat("found", length(seus), "hypoxia objects\n")
seus <- head(seus, 3)
cat("smoke samples:\n"); cat(paste0("  ", seus, collapse = "\n"), "\n")

# --- 1. does the object carry what the replay needs? ------------------------
s <- readRDS(seus[[1]])
cat("\n--- object check:", basename(seus[[1]]), "---\n")
cat("cells:      ", ncol(s), "\n")
cat("assays:     ", paste(names(s@assays), collapse = ", "), "\n")
cat("reductions: ", paste(names(s@reductions), collapse = ", "), "\n")
cat("score cols: ", paste(intersect(c("hypoxia", "MT", "hypoxia_score", "Phase"),
                                    colnames(s@meta.data)), collapse = ", "), "\n")
rm(s); invisible(gc())

# --- 2. replay + reconcile + render -----------------------------------------
out_pdf <- "results/hypoxia_cluster_split/SMOKE_hypoxia_mito_score_boxplots.pdf"
cat("\n--- plot_hypoxia_mito_score_boxplots ---\n")
res <- plot_hypoxia_mito_score_boxplots(
  split_log_csv     = log_csv,
  hypoxia_seu_paths = seus,
  out_pdf           = out_pdf,
  resolutions       = seq(0.2, 1.2, by = 0.2),
  split_assay       = "gene")

cat("\nreturned: ", res, "\n")
cat("pdf ok:   ", !is.na(res) && file.exists(res), "\n")
if (!is.na(res) && file.exists(res))
  cat("pdf size: ", round(file.size(res) / 1024), " KB\n")

# --- 3. show the reconciliation explicitly ----------------------------------
# Re-run the join here so the numbers are visible rather than only messaged.
cat("\n--- explicit reconciliation ---\n")
d <- readr::read_csv(log_csv, show_col_types = FALSE) |>
  mutate(resolution = round(as.numeric(resolution), 3)) |>
  filter(round == 1L, sample_id %in% stringr::str_extract(seus, "SR[RX][0-9]+")) |>
  mutate(cluster = as.character(cluster))

r <- replay_cluster_score_means(seus, resolutions = seq(0.2, 1.2, by = 0.2),
                                split_assay = "gene") |>
  mutate(cluster = as.character(cluster))

j <- inner_join(d, r, by = c("sample_id", "resolution", "cluster"),
                suffix = c("", "_replay"))
cat("logged round-1 clusters for these samples:", nrow(d), "\n")
cat("joined:                                   ", nrow(j), "\n")
cat("n_cells identical:                        ", sum(j$n_cells == j$n_cells_replay), "\n")
cat("mean_score within 1e-6:                   ",
    sum(abs(j$mean_score - j$mean_score_replay) < 1e-6), "\n")
cat("max |mean_score diff|:                    ",
    format(max(abs(j$mean_score - j$mean_score_replay)), digits = 3), "\n")

bad <- j[j$n_cells != j$n_cells_replay |
           abs(j$mean_score - j$mean_score_replay) >= 1e-6, ]
if (nrow(bad) > 0) {
  cat("\n!! ", nrow(bad), " clusters did NOT reconcile:\n")
  print(head(bad[, c("sample_id", "resolution", "cluster", "n_cells",
                     "n_cells_replay", "mean_score", "mean_score_replay")], 10))
} else {
  cat("\nall joined clusters reconcile exactly\n")
}

cat("\nSMOKE DONE\n")
