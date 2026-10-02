#!/usr/bin/env Rscript
# Build numbat inputs for one tumor from its technical-replicate captures, in the
# formats pipeline/scripts/run_numbat.R already reads (so it needs no changes):
#   <outdir>/matrix/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}  -- Read10X dir
#   <outdir>/seu_meta.rds     -- minimal Seurat object: just the metadata
#                                (rownames + `type`) that run_numbat.R filters on
#   <outdir>/allele_counts.tsv.gz
#
# The captures are independent 10x runs (~1% barcode overlap, i.e. whitelist
# collisions), so cells are concatenated, never summed, and every barcode is
# prefixed "<SRX>_" in all three files.
#
# Expects the jointly phased per-replicate tables from
# src/phase_merged_replicates.R at <outdir>/<SRX>_allele_counts.tsv.gz.
#
# Usage: Rscript src/build_merged_replicate_inputs.R <SRX_a,SRX_b> <outdir>

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
samples <- strsplit(args[[1]], ",")[[1]]
outdir <- args[[2]]
proj <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj"
mtx_dir <- function(s) file.path(proj, "output/cellranger", s, "outs/filtered_feature_bc_matrix")

## expression: cbind raw 10x matrices with prefixed barcodes
features <- lapply(samples, function(s) fread(file.path(mtx_dir(s), "features.tsv.gz"), header = FALSE))
stopifnot(all(vapply(features[-1], identical, logical(1), features[[1]])))

mats <- lapply(samples, function(s) {
  m <- readMM(file.path(mtx_dir(s), "matrix.mtx.gz"))
  colnames(m) <- paste0(s, "_", readLines(file.path(mtx_dir(s), "barcodes.tsv.gz")))
  as(m, "CsparseMatrix")
})
count_mat <- do.call(cbind, mats)
stopifnot(!anyDuplicated(colnames(count_mat)))

dir.create(file.path(outdir, "matrix"), recursive = TRUE, showWarnings = FALSE)
writeMM(count_mat, file.path(outdir, "matrix/matrix.mtx"))
system2("gzip", c("-f", file.path(outdir, "matrix/matrix.mtx")))
gz <- gzfile(file.path(outdir, "matrix/barcodes.tsv.gz"), "w")
writeLines(colnames(count_mat), gz)
close(gz)
fwrite(features[[1]], file.path(outdir, "matrix/features.tsv.gz"), sep = "\t", col.names = FALSE)

## metadata: keep each replicate's own cell filter + `type`, prefixed
meta <- rbindlist(lapply(samples, function(s) {
  md <- readRDS(file.path(proj, "output/seurat", paste0(s, "_seu.rds")))@meta.data
  cells <- sub("\\.(\\d+)$", "-\\1", rownames(md))
  data.table(
    cell = paste0(s, "_", cells),
    replicate = s,
    type = if ("type" %in% names(md)) as.character(md$type) else NA_character_
  )
}))
stopifnot(!anyDuplicated(meta$cell))
meta_df <- data.frame(meta[, .(replicate, type)], row.names = meta$cell)
dummy <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(1, nrow(meta_df)),
                              dimnames = list("placeholder", meta$cell))
seu <- SeuratObject::CreateSeuratObject(counts = dummy, meta.data = meta_df)
stopifnot(identical(rownames(seu@meta.data), meta$cell))
saveRDS(seu, file.path(outdir, "seu_meta.rds"))

## alleles: prefix and stack the jointly phased tables
allele <- rbindlist(lapply(samples, function(s) {
  a <- fread(file.path(outdir, paste0(s, "_allele_counts.tsv.gz")))
  a[, cell := paste0(s, "_", cell)]
}))
stopifnot(all(unique(allele$cell) %in% colnames(count_mat)))
fwrite(allele, file.path(outdir, "allele_counts.tsv.gz"), sep = "\t")

cat(sprintf("cells: matrix=%d  seu=%d  allele=%d | het SNPs=%d\n",
            ncol(count_mat), nrow(meta), uniqueN(allele$cell), uniqueN(allele$snp_id)))
print(meta[, .N, by = replicate])
