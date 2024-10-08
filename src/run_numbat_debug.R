#!/usr/bin Rscript

if(!exists("expression_only")){
  expression_only <- FALSE
}

# Rprof(rprof_out)

print(paste0("num_clones: ", init_k))

library(numbat)
# undebug(run_numbat)
library(Seurat)
library(readr)
library(magrittr)
conflicted::conflict_prefer("rowSums", "Matrix")

out_dir = "output/numbat_sridhar/SRR13884242/"

matrix_file = "output/cellranger/SRR13884242/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
print(matrix_file)

matrix_dir = fs::path_dir(matrix_file)

count_mat <- Seurat::Read10X(matrix_dir)

# subset count matrix by cell to ease tree construction ------------------------------
cell_ceiling = 1e5
cell_ceiling <- as.numeric(cell_ceiling)
if(cell_ceiling > 0){

	cell_ceiling = ifelse(ncol(count_mat) > cell_ceiling, cell_ceiling, ncol(count_mat))

	count_mat <- count_mat[,sample(ncol(count_mat),size=cell_ceiling,replace=FALSE)]
}

# test dropping read numbers ------------------------------
if(!as.numeric(read_prop) == 1){
	count_mat <- scuttle::downsampleMatrix(count_mat, as.numeric(read_prop))
}
allele_df = "output/numbat/SRR13884242/SRR13884242_allele_counts.tsv.gz"
allele_df = data.table::fread(allele_df) %>%
    tidyr::drop_na(CHROM) %>% # drop X chromosome values
  identity()

seu_path = "output/seurat/SRR13884242_seu.rds"
myseu <- readRDS(seu_path)

count_mat <- count_mat[,colnames(count_mat) %in% colnames(myseu)]

ref_path = "~/Homo_sapiens/numbat/sridhar_ref.rds"
print(ref_path)
ref_internal <- readRDS(ref_path)

bulk = numbat:::get_bulk(
	count_mat = count_mat,
	lambdas_ref = ref_internal,
	df_allele = allele_df,
	gtf = gtf_hg38
)

segs_loh = bulk %>% numbat:::detect_clonal_loh(t = as.numeric(1e-05))

# segs_loh = NULL

# run
out = numbat::run_numbat(
	count_mat, # gene x cell integer UMI count matrix done
	lambdas_ref = ref_internal, # reference expression profile, a gene x cell type normalized expression level matrix done
	allele_df, # allele dataframe generated by pileup_and_phase script done
	genome = "hg38",
	min_cells = 10,
	segs_loh = segs_loh,
	multi_allelic = FALSE,
	out_dir = out_dir,
	expression_only = FALSE,
	t = 1e-05,
	alpha = 1e-04,
	gamma = 50,
	init_k = 5,
	max_cost = 2425.5,
	n_cut = 0,
	max_iter = 1,
	max_nni = 100,
	min_depth = 0,
	use_loh = "auto",
	min_LLR = 2,
	min_overlap = 0.45,
	max_entropy = 0.8,
	skip_nj = TRUE,
	diploid_chroms = NULL,
	ncores = 1,
	ncores_nni = 1,
	common_diploid = TRUE,
	tau = 0.3,
	check_convergence = FALSE,
	plot = TRUE
	)

# settings ------------------------------
out_dir = "output/numbat_sridhar/SRR13884242/"
i=2
exp_post <- fread(glue("{out_dir}/exp_post_{i}.tsv"),
       sep = "\t")
allele_post <- fread(glue("{out_dir}/allele_post_{i}.tsv"),
       sep = "\t")
joint_post <- fread(glue("{out_dir}/joint_post_{i}.tsv"),
       sep = "\t")
tree_list <- readRDS(glue("{out_dir}/tree_list_{i}.rds"))
treeML <- readRDS(glue("{out_dir}/treeML_{i}.rds"))
gtree <- readRDS(glue("{out_dir}/tree_final_{i}.rds"))
G_m <- readRDS(glue("{out_dir}/mut_graph_{i}.rds"))
clone_post = fread("output/numbat_sridhar/SRR13884242/clone_post_2.tsv", sep = "\t")
subtrees = fread("output/numbat_sridhar/SRR13884242/bulk_subtrees_2.tsv.gz", sep = "\t")
clones <- readRDS("output/numbat_sridhar/SRR13884242/clones_2.rds")
segs_consensus = fread("output/numbat_sridhar/SRR13884242/segs_consensus_2.tsv", sep = "\t")

done_file <- fs::path(out_dir, "done.txt")

write(ref_path, file = done_file)
