#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

tar_load("large_clone_simplifications")

large_clone_simplifications %>%
  map(tibble::enframe, "scna", "seg") %>%
  dplyr::bind_rows(.id= "sample_id") %>%
  write_csv("results/large_clone_simplifications.csv")

new_adata_filtered_meta <- read_csv("results/new_adata_filtered_meta.csv")

large_clone_simplifications  = read_csv("results/large_clone_simplifications.csv")

new_meta <-
  new_adata_filtered_meta %>%
  dplyr::mutate(seg = str_split(GT_opt, pattern = ",")) %>%
  # select(seg) %>%
  tidyr::unnest(seg) %>%
  dplyr::left_join(large_clone_simplifications, by = c("sample_id", "seg")) %>%
  # dplyr::select(-seg) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::group_by(sample_id, GT_opt) %>%
  dplyr::distinct(sample_id, GT_opt, scna) %>%
  dplyr::filter(!is.na(scna)) %>%
  dplyr::summarize(scna = toString(scna)) %>%
  identity()

new_meta <-
  new_adata_filtered_meta %>%
  dplyr::left_join(new_meta, by = c("sample_id", "GT_opt")) %>%
  identity()
