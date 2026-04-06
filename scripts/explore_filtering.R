#!/usr/bin/env Rscript
# Interactive script for exploring filtering threshold effects.
#
# Usage (from project root, interactively in R):
#   source("scripts/explore_filtering.R")
#
# Or run a specific section by copying the desired block into an R session
# after tar_load(filter_inspection_metadata).

library(targets)

# Load the per-cell metadata tibble cached by the pipeline.
# This only reads targets store metadata — no Seurat RDS is loaded here.
tar_load(filter_inspection_metadata, store = "_targets_r431")
meta <- dplyr::bind_rows(filter_inspection_metadata)

message("Loaded ", nrow(meta), " cells across ", dplyr::n_distinct(meta$sample_id), " samples.")
message("Samples: ", paste(unique(meta$sample_id), collapse = ", "))

# ── Single-sample threshold sweep ────────────────────────────────────────────
# Visualise how cell retention varies across QC threshold combinations for one sample.
p_sweep <- plot_filter_sweep(dplyr::filter(meta, sample_id == "SRR14800534"))
print(p_sweep)

# ── Custom threshold grid ─────────────────────────────────────────────────────
# Use a narrower grid when you have a hypothesis about the relevant range.
custom_grid <- expand.grid(
  mito     = c(5, 10, 20),
  nCount   = 1000,
  nFeature = 1000,
  stringsAsFactors = FALSE
)
p_custom <- plot_filter_sweep(meta, param_grid = custom_grid)
print(p_custom)

# ── Staged summary at a specific threshold set ────────────────────────────────
# Inspect how many cells are removed at each filtering stage.
staged <- apply_filter_criteria(
  meta,
  mito_threshold     = 15,
  nCount_threshold   = 500,
  nFeature_threshold = 500
)
print(staged)

# ── All-sample sweep (combined) ───────────────────────────────────────────────
# Overlay all samples in one plot (colour = scna, facet = qc_metric).
p_all <- plot_filter_sweep(meta)
print(p_all)
