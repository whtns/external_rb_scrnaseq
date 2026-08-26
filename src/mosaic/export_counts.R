#!/usr/bin/env Rscript
# Issue #40, Stage B. Export raw integer counts (+ cell metadata) for ONE
# paper_retained sample, for downstream mosaicMPI cNMF. Writes a sparse MatrixMarket
# triple (features x cells) plus barcodes/genes/metadata, so build_h5ad.py can
# construct the mosaicMPI Dataset without a multi-GB dense-TSV round-trip.
#
# Correctness (per the phone-a-friend review): cNMF needs RAW COUNTS. We pull the
# explicit counts layer, assert it is integer-valued, and confirm rownames are gene
# SYMBOLS (the reference GMT is symbol-keyed). No normalization here.
#
# Usage: Rscript src/mosaic/export_counts.R <SRX> [out_dir]
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: export_counts.R <SRX> [out_dir]")
srx     <- args[[1]]
out_dir <- if (length(args) >= 2) args[[2]] else file.path("output/mosaicmpi_rb", srx)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seu_path <- file.path("output/seurat", paste0(srx, "_filtered_seu.rds"))
if (!file.exists(seu_path)) stop("missing Seurat object: ", seu_path)
seu <- readRDS(seu_path)

# Prefer the gene assay (full raw gene set); fall back to RNA. NOT SCT (subset +
# corrected, not raw counts).
assay <- intersect(c("gene", "RNA"), names(seu@assays))
if (length(assay) == 0) stop("no gene/RNA assay in ", seu_path,
                             " (have: ", paste(names(seu@assays), collapse=", "), ")")
assay <- assay[[1]]
counts <- SeuratObject::GetAssayData(seu, assay = assay, layer = "counts")
if (is.null(counts) || length(counts@x) == 0)
  stop("empty counts layer for assay ", assay, " in ", seu_path)

# Assertions: integer counts, symbol rownames.
if (!all(counts@x == floor(counts@x)))
  stop("counts for ", srx, " assay ", assay, " are not integer-valued -- refusing to ",
       "export normalized/corrected data as cNMF input")
rn <- rownames(counts)
looks_ensembl <- mean(grepl("^ENSG[0-9]{11}", rn)) > 0.5
if (looks_ensembl)
  stop("rownames of ", srx, " look like Ensembl IDs, not symbols; the reference GMT is ",
       "symbol-keyed. Map to symbols before export.")

# Metadata layers useful downstream (per-cell score comparison + grouping).
keep_meta <- intersect(
  c("hypoxia_score", "S.Score", "G2M.Score", "Phase", "clusters", "clone_opt",
    "SCT_snn_res.0.6"),
  colnames(seu@meta.data))
meta <- seu@meta.data[, keep_meta, drop = FALSE]
meta <- data.frame(barcode = rownames(meta), meta, check.names = FALSE,
                   row.names = NULL)

Matrix::writeMM(counts, file.path(out_dir, "counts.mtx"))
writeLines(colnames(counts), file.path(out_dir, "barcodes.txt"))
writeLines(rn,               file.path(out_dir, "genes.txt"))
readr::write_tsv(meta,       file.path(out_dir, "metadata.tsv"))

cat(sprintf("exported %s: assay=%s, %d genes x %d cells, integer counts, symbol rownames\n",
            srx, assay, nrow(counts), ncol(counts)))
cat("  ->", out_dir, "(counts.mtx, barcodes.txt, genes.txt, metadata.tsv)\n")
