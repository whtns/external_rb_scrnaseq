#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
library(patchwork)
# asdf

bulk_subtrees_2_4242 <- read_tsv("output/numbat_sridhar/SRR13884242/bulk_subtrees_2.tsv.gz") %>%
  dplyr::mutate(CHROM = factor(CHROM))

bulk_subtrees_2_4242_retest_1 <- read_tsv("output/numbat_sridhar/SRR13884242/bulk_subtrees_retest_1.tsv.gz") %>%
  dplyr::mutate(CHROM = factor(CHROM))

bulk_subtrees_2_4242_retest_2 <- read_tsv("output/numbat_sridhar/SRR13884242/bulk_subtrees_retest_2.tsv.gz") %>%
  dplyr::mutate(CHROM = factor(CHROM))

bulk_clones_final_4242 <- read_tsv("output/numbat_sridhar/SRR13884242/bulk_clones_final.tsv.gz") %>%
  dplyr::mutate(CHROM = factor(CHROM))

# debug(numbat:::run_group_hmms)

bulk_subtrees <- numbat:::run_group_hmms(bulk_subtrees_2_4242)

bulk_subtrees_r1 <- numbat:::run_group_hmms(bulk_subtrees_2_4242_retest_1)

bulk_subtrees_r2 <- numbat:::run_group_hmms(bulk_subtrees_2_4242_retest_2)

bulk_clones <- numbat:::run_group_hmms(bulk_clones_final_4242)

p_subtrees_r1 = numbat:::plot_bulks(bulk_subtrees_r1, min_LLR = 2, use_pos = TRUE, genome = 'hg38')

pdf("results/SRR13884242_clone_num_debug_plots.pdf", width = 8, height = 6)
p_subtrees = numbat:::plot_bulks(bulk_subtrees, min_LLR = 2, use_pos = TRUE, genome = 'hg38') +
  plot_annotation(title = "subtrees 2")
p_subtrees

p_subtrees_r2 = numbat:::plot_bulks(bulk_subtrees_r2, min_LLR = 2, use_pos = TRUE, genome = 'hg38') +
  plot_annotation(title = "subtrees 2 retest")
p_subtrees_r2

p_clones_final = numbat:::plot_bulks(bulk_clones, min_LLR = 2, use_pos = TRUE, genome = 'hg38') +
  plot_annotation(title = "clones final")
p_clones_final

dev.off()

# test run_numbat ------------------------------

# run_numbat()
