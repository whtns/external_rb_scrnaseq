#!/usr/bin/env Rscript
# Issue #40 (interpretation add-on). Build a BROAD reference GMT: all 50 MSigDB
# HALLMARK sets. Unlike data/mosaic_reference_programs.gmt (our 4 hand-picked
# hypoxia/cell-cycle signatures -- the confirmatory view), this lets each mosaicMPI
# community be annotated data-driven, against a library we did NOT curate for this
# question. This is what makes the large pan-tumor community (#1), which our 4
# references cannot see, interpretable.
#
# GMT line: <set_name>\t<description>\t<gene1>\t<gene2>...  (symbols only)
# Usage: Rscript src/mosaic/build_hallmark_gmt.R [out_gmt]
suppressPackageStartupMessages({ library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
out_gmt <- if (length(args) >= 1) args[[1]] else "data/mosaic_hallmark.gmt"

h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
sets <- split(h$gene_symbol, h$gs_name)
sets <- lapply(sets, function(g) unique(g[nzchar(g) & !is.na(g)]))
sets <- sets[vapply(sets, length, integer(1)) >= 5]           # ssGSEA min set size
stopifnot(length(sets) >= 45)                                  # expect ~50 Hallmark sets

dir.create(dirname(out_gmt), showWarnings = FALSE, recursive = TRUE)
lines <- vapply(names(sets), function(nm)
  paste(c(nm, "hallmark", sets[[nm]]), collapse = "\t"), character(1))
writeLines(lines, out_gmt)
cat("wrote", out_gmt, "with", length(sets), "HALLMARK sets\n")
