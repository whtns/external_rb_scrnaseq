#!/usr/bin/env Rscript
# Single-branch effect_of_filtering with the new cluster-based hypoxia split.
# Mirrors the effect_of_filtering target body for one sample, using the
# *_hypoxia_low/high_seu.rds produced by split_hypoxia_by_clusters().
suppressPackageStartupMessages(source("packages.R"))
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
sample_id <- if (length(args) >= 1) args[[1]] else "SRX10031194"
store <- "_targets_r431"

unfiltered_path   <- glue::glue("output/seurat/{sample_id}_unfiltered_seu.rds")
filtered_path     <- glue::glue("output/seurat/{sample_id}_filtered_seu.rds")
low_hypoxia_path  <- glue::glue("output/seurat/{sample_id}_hypoxia_low_seu.rds")
high_hypoxia_path <- glue::glue("output/seurat/{sample_id}_hypoxia_high_seu.rds")

stopifnot(file.exists(unfiltered_path), file.exists(filtered_path),
          file.exists(low_hypoxia_path), file.exists(high_hypoxia_path))

# Reuse the dictionary already computed by the pipeline (no rebuild).
cluster_dictionary <- targets::tar_read(cluster_dictionary, store = store)

cat("== running plot_effect_of_filtering for", sample_id, "==\n")
out <- plot_effect_of_filtering(
  unfiltered_path, filtered_path,
  cluster_dictionary    = cluster_dictionary,
  low_hypoxia_seu_path  = low_hypoxia_path,
  high_hypoxia_seu_path = high_hypoxia_path
)
cat("\n== wrote:", out, "==\n")
cat("exists:", file.exists(out), " size:",
    if (file.exists(out)) file.size(out) else NA, "\n")
