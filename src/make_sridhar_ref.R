#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)

sridhar_seu <- readRDS("/project2/cobrinik_1090/external_rb_scrnaseq_proj/data/sridhar_unfiltered_seu.rds")

refs <- list(
plae_ref = readRDS("/project2/cobrinik_1090/Homo_sapiens/numbat/plae_ref.rds"), # plae input
sridhar_ref = readRDS("/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds"), # based off data
new_sridhar_ref = readRDS("/project2/cobrinik_1090/Homo_sapiens/numbat/new_sridhar_ref.rds") # based off counts
)

nb <- readRDS("output/numbat_sridhar/SRX11133594_numbat.rds")

nb$plot_exp_roll(
  k = 5,
  n_sample = 300,  # Explicitly match resolution
  lim = 0.8         # Hard-cap the color contrast scale (lower = more vibrant)
)

fetal_clustifyr <- readRDS("/project2/cobrinik_1090/Homo_sapiens/numbat/fetal_clustifyr_seu.rds")



# count_mat is a gene x cell raw count matrices
count_mat <- Seurat::GetAssayData(fetal_clustifyr, layer = "counts", assay = "gene")
# cell_annot is a dataframe with columns "cell" and "group"
cell_annot <- fetal_clustifyr@meta.data["clustifyr_type"] %>%
	# as.data.frame() %>%
	tibble::rownames_to_column("cell") %>%
	dplyr::rename(group = clustifyr_type) %>%
	dplyr::mutate(group = ifelse(is.na(group), "NA", group))

ref_internal = aggregate_counts(count_mat, cell_annot)

ref_internal <- ref_internal[,colnames(ref_internal) != "NA"]

# saveRDS(ref_internal, "~/Homo_sapiens/numbat/sridhar_ref.rds")
# saveRDS(ref_internal, "/project2/cobrinik_1090/Homo_sapients/new_sridhar_ref.rds")

new_sridhar <- readRDS("/project2/cobrinik_1090/Homo_sapiens/numbat/new_sridhar_ref.rds")

canonical_sridhar <- readRDS("/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds")
