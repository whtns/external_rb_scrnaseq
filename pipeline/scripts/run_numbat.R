#!/usr/bin/env Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

if(!exists("expression_only")){
  expression_only <- FALSE
}

if(!exists("subset_bad_cell_types")){
  subset_bad_cell_types <- TRUE
}

if(!exists("min_allele_depth")){
	min_allele_depth <- 5
}

if(!exists("min_snps_per_cell")){
	min_snps_per_cell <- 50
}

if(!exists("max_nni")){
	max_nni <- 100
}

if(!exists("retry_without_nni")){
	retry_without_nni <- TRUE
}

parse_bool <- function(x, default = TRUE) {
	if (is.logical(x)) return(x)
	if (is.numeric(x)) return(x != 0)
	if (is.null(x)) return(default)
	x <- tolower(trimws(as.character(x)))
	if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
	if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
	default
}

parse_num <- function(x, default) {
	if (is.null(x)) return(default)
	if (is.numeric(x)) return(as.numeric(x))
	v <- suppressWarnings(as.numeric(as.character(x)))
	if (is.na(v)) return(default)
	v
}

subset_bad_cell_types <- parse_bool(subset_bad_cell_types, TRUE)
retry_without_nni <- parse_bool(retry_without_nni, TRUE)
min_allele_depth <- parse_num(min_allele_depth, 5)
min_snps_per_cell <- parse_num(min_snps_per_cell, 50)
max_nni <- as.integer(parse_num(max_nni, 100))

# Rprof(rprof_out)

print(paste0("num_clones: ", init_k))

library(numbat)
library(Seurat)
library(readr)
library(magrittr)
conflicted::conflict_prefer("rowSums", "Matrix")
bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))

print(matrix_file)

matrix_dir = fs::path_dir(matrix_file)

count_mat <- Seurat::Read10X(matrix_dir)

myseu <- readRDS(seu_path)

seu_cells <- rownames(myseu@meta.data)
if (is.null(seu_cells) || length(seu_cells) == 0) {
	stop("No cell IDs found in Seurat metadata rownames.")
}

# Filter known non-tumor cell types only when requested and available.
if (subset_bad_cell_types && "type" %in% colnames(myseu@meta.data)) {
	cell_type <- myseu@meta.data[["type"]]
	names(cell_type) <- seu_cells
	keep_cells <- !(cell_type %in% bad_cell_types)
	keep_cells[is.na(keep_cells)] <- TRUE
	seu_cells <- names(keep_cells)[keep_cells]
	if (length(seu_cells) == 0) {
		warning("Type-based filtering removed all cells; falling back to all metadata cells.")
		seu_cells <- rownames(myseu@meta.data)
	}
} else if (!subset_bad_cell_types) {
	warning("Skipping bad_cell_types filtering because subset_bad_cell_types is FALSE.")
} else {
	warning("Metadata column 'type' not found; skipping bad_cell_types filtering.")
}

# drop cells absent from metadata-filtered set
count_mat <- count_mat[,colnames(count_mat) %in% seu_cells]

# subset count matrix by cell to ease tree construction ------------------------------
# cell_ceiling <- as.numeric(cell_ceiling)
# if(cell_ceiling > 0){
# 
# 	cell_ceiling = ifelse(ncol(count_mat) > cell_ceiling, cell_ceiling, ncol(count_mat))
# 
# 	count_mat <- count_mat[,sample(ncol(count_mat),size=cell_ceiling,replace=FALSE)]
# }

# test dropping read numbers ------------------------------
if(!as.numeric(read_prop) == 1){
	count_mat <- scuttle::downsampleMatrix(count_mat, as.numeric(read_prop))
}

allele_df = data.table::fread(allele_df) |> 
	tidyr::drop_na(CHROM) |> # drop X chromosome values
	dplyr::filter(cell %in% seu_cells) |> # drop cells absent from metadata-filtered set
	identity()

if (!"cell" %in% colnames(allele_df)) {
	stop("Allele table does not contain required column 'cell'.")
}

if (ncol(count_mat) == 0 || nrow(allele_df) == 0) {
	stop("Empty input after initial filtering: count_mat columns=", ncol(count_mat), ", allele rows=", nrow(allele_df))
}

shared_cells <- intersect(colnames(count_mat), unique(allele_df$cell))
if (length(shared_cells) == 0) {
	stop("No overlapping cell barcodes between expression and allele data after filtering.")
}

if (length(shared_cells) < 100) {
	warning("Low overlap between expression and allele data (", length(shared_cells), " cells).")
}

count_mat <- count_mat[, colnames(count_mat) %in% shared_cells, drop = FALSE]
allele_df <- dplyr::filter(allele_df, cell %in% shared_cells)

ref_col <- intersect(c("n_ref", "ref_count", "REF_COUNT"), colnames(allele_df))
alt_col <- intersect(c("n_alt", "alt_count", "ALT_COUNT"), colnames(allele_df))

if (length(ref_col) > 0 && length(alt_col) > 0) {
	ref_col <- ref_col[[1]]
	alt_col <- alt_col[[1]]

	depth <- suppressWarnings(as.numeric(allele_df[[ref_col]]) + as.numeric(allele_df[[alt_col]]))
	keep_depth <- is.finite(depth) & depth >= min_allele_depth
	allele_df <- allele_df[keep_depth, , drop = FALSE]

	if (nrow(allele_df) == 0) {
		stop("All allele rows removed by min_allele_depth filter (", min_allele_depth, ").")
	}
} else {
	warning("Could not find n_ref/n_alt-style columns; skipping allele depth filter.")
}

cell_snp_counts <- table(allele_df$cell)
cells_with_snps <- names(cell_snp_counts[cell_snp_counts >= min_snps_per_cell])
if (length(cells_with_snps) == 0) {
	warning("No cells met min_snps_per_cell=", min_snps_per_cell, "; keeping all overlapping cells.")
	cells_with_snps <- shared_cells
}

count_mat <- count_mat[, colnames(count_mat) %in% cells_with_snps, drop = FALSE]
allele_df <- dplyr::filter(allele_df, cell %in% cells_with_snps)

if (ncol(count_mat) == 0 || nrow(allele_df) == 0) {
	stop("Empty input after SNP filtering: count_mat columns=", ncol(count_mat), ", allele rows=", nrow(allele_df))
}

if (anyNA(allele_df$cell) || anyNA(allele_df$CHROM)) {
	stop("Unexpected NA values in required allele columns after filtering.")
}

print(paste0("Cells after overlap/SNP filters: ", ncol(count_mat)))
print(paste0("Allele rows after filters: ", nrow(allele_df)))

# assemble a normal reference from non-RB cell types ------------------------------

# num_out_cells <- sum(myseu$type %in% bad_cell_types)

# if(num_out_cells > 100){
# 	normal_seu = myseu[,myseu$type %in% bad_cell_types]
#
# 	normal_reference_mat <-
# 		normal_seu %>%
# 		GetAssayData(slot = "counts")
#
# 	normal_cell_annot <-
# 		normal_seu@meta.data[c("type")] %>%
# 		tibble::rownames_to_column("cell") %>%
# 		dplyr::rename(group = type) %>%
# 		identity()
#
# 	ref_internal = numbat::aggregate_counts(normal_reference_mat, normal_cell_annot)
# } else {
# 	ref_internal <- readRDS(ref_path)
# }

print(ref_path)
ref_internal <- readRDS(ref_path)

bulk = numbat:::get_bulk(
	count_mat = count_mat,
	lambdas_ref = ref_internal,
	df_allele = allele_df,
	gtf = gtf_hg38
)

segs_loh = bulk %>% numbat:::detect_clonal_loh(t = as.numeric(t))

# segs_loh = NULL

# run
run_numbat_once <- function(max_nni_local) {
	numbat::run_numbat(
		count_mat, # gene x cell integer UMI count matrix done
		lambdas_ref = ref_internal, # reference expression profile, a gene x cell type normalized expression level matrix done
		allele_df, # allele dataframe generated by pileup_and_phase script done
		genome = "hg38",
		min_cells = 10,
		t = as.numeric(t),
		gamma = as.numeric(gamma),
		max_iter = as.integer(max_iter),
		max_entropy = max_entropy,
		segs_loh = segs_loh,
		min_LLR = min_LLR,
		init_k = init_k,
		ncores = as.integer(ncores),
		plot = TRUE,
		multi_allelic = FALSE,
		common_diploid = TRUE,
		out_dir = out_dir,
		tau = as.numeric(tau),
		skip_nj = TRUE,
		max_nni = as.integer(max_nni_local)
	)
}

out <- tryCatch(
	run_numbat_once(max_nni),
	error = function(e) {
		msg <- conditionMessage(e)
		if (retry_without_nni && grepl("perform_nni|missing value where TRUE/FALSE needed|arguments imply differing number of rows", msg, ignore.case = TRUE)) {
			warning("run_numbat failed during NNI optimization; retrying with max_nni=0.")
			return(run_numbat_once(0L))
		}
		stop(e)
	}
)

done_file <- fs::path(out_dir, "done.txt")

write(ref_path, file = done_file)
