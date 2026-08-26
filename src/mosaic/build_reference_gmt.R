#!/usr/bin/env Rscript
# Issue #40, Stage A. Build the reference-program GMT that mosaicMPI programs get
# scored against via `mosaicmpi ssgsea`. These are the "programs we already use":
#   HALLMARK_HYPOXIA   - the hypoxia signature (msigdbr; heatmap_functions.R)
#   TIROSH_S_PHASE     - Tirosh S-phase cell-cycle markers (Seurat::cc.genes.updated.2019)
#   TIROSH_G2M         - Tirosh G2/M markers
#   GIOTTI_CELL_CYCLE  - the 16 giotti CC categories collapsed to ONE set (headline)
#   GIOTTI__<term>     - the 16 categories kept individually (for the CSV / detail)
#
# GMT line: <set_name>\t<description>\t<gene1>\t<gene2>...
# Symbols only -- mosaicMPI programs use gene-symbol var_names.
#
# Usage: Rscript src/mosaic/build_reference_gmt.R [out_gmt]
suppressPackageStartupMessages({
  library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
out_gmt <- if (length(args) >= 1) args[[1]] else "data/mosaic_reference_programs.gmt"

sets <- list()

## HALLMARK_HYPOXIA -- exact set the pipeline uses (heatmap_functions.R).
hy <- msigdbr::msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_name == "HALLMARK_HYPOXIA") %>%
  dplyr::pull(gene_symbol) %>% unique()
stopifnot(length(hy) > 100)  # HALLMARK_HYPOXIA is ~200 genes; guard against an empty pull
sets[["HALLMARK_HYPOXIA"]] <- hy

## Tirosh S and G2M markers, as used by Seurat::CellCycleScoring throughout the paper.
cc <- Seurat::cc.genes.updated.2019
sets[["TIROSH_S_PHASE"]] <- unique(cc$s.genes)
sets[["TIROSH_G2M"]]     <- unique(cc$g2m.genes)

## Giotti cell-cycle categories (data/giotti_cc_genes.tsv). One combined set for the
## headline heatmap (the 16 categories are overlapping CC subsets), plus each
## category individually for the detail table.
gio <- readr::read_tsv("data/giotti_cc_genes.tsv", show_col_types = FALSE)
sets[["GIOTTI_CELL_CYCLE"]] <- unique(gio$symbol)
for (tm in sort(unique(gio$term))) {
  key <- paste0("GIOTTI__", gsub("[^A-Za-z0-9]+", "_", toupper(tm)))
  sets[[key]] <- unique(gio$symbol[gio$term == tm])
}

## Write GMT.
lines <- vapply(names(sets), function(nm) {
  paste(c(nm, nm, sets[[nm]]), collapse = "\t")
}, character(1))
dir.create(dirname(out_gmt), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, out_gmt)

cat("wrote", out_gmt, "with", length(sets), "gene sets:\n")
for (nm in names(sets)) cat(sprintf("  %-40s %d genes\n", nm, length(sets[[nm]])))
