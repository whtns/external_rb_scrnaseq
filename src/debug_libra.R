#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(Libra)
library(Seurat)

inspect_libra_output <- function(df){

  select_annotables <-
    annotables::grch38 %>%
    dplyr::select(symbol, description) %>%
    dplyr::distinct(.keep_all = TRUE)

  df %>%
    dplyr::filter(p_val < 0.05) %>%
    dplyr::group_by(cell_type) %>%
    arrange(desc(abs(avg_logFC))) %>%
    # slice_head(n=20) %>%
    dplyr::left_join(select_annotables, by = c("gene" = "symbol"))
}


seu <- readRDS("output/seurat/CopyOfSRR14800534_SRR14800535_SRR14800536_seu.rds")
seu <- seu[,!is.na(seu$clone_opt)]

dropped_cell_types <- c("APOE", "MALAT1", "rod")
seu <- seu[,(!seu$abbreviation %in% dropped_cell_types)]

seu$abbreviation <- case_when(seu$abbreviation %in% c("ARL1IP1", "HIST1H4C") ~ "G2M",
          seu$abbreviation %in% c("PCLAF") ~ "TFF1",
          TRUE ~ seu$abbreviation)

# seu$abbreviation[!seu$abbreviation %in% c("G2M", "TFF1")] <- "early_RB"
seu@meta.data["cell_type"] = factor(seu@meta.data$abbreviation)

seu$clone_opt[seu$clone_opt == "4"] <- "3"
seu@meta.data["label"] = factor(seu@meta.data$clone_opt)

seu@meta.data["replicate"] = factor(seu@meta.data$batch)

expr_input = GetAssayData(seu, assay = "gene", slot = "counts") %>%
  # na.omit() %>%
  identity()

# single cell
singlecell_DE = run_de(expr_input,
                       meta = seu@meta.data,
                       de_family = "singlecell",
                       de_method = "wilcox",
                       min_reps = 1)

inspect_singlecell_DE <-
singlecell_DE %>%
  inspect_libra_output()

# pseudobulk 1 v. 2 ------------------------------

seu_12 <- seu[,seu$clone_opt %in% c("1", "2")]

expr_input = GetAssayData(seu_12, assay = "gene", slot = "counts") %>%
  # na.omit() %>%
  identity()

pseudobulk_DE_12 = run_de(expr_input,
                       meta = seu_12@meta.data,
                       de_family = "pseudobulk",
                       de_method = "DESeq2",
                       de_type = "LRT",
                       min_reps = 2,
                       min_cells = 2)

inspect_pseudobulk_DE_12 <-
  pseudobulk_DE_12 %>%
  inspect_libra_output()

# pseudobulk 2 v. 3 ------------------------------

seu_23 <- seu[,seu$clone_opt %in% c("2", "3")]

expr_input = GetAssayData(seu_23, assay = "gene", slot = "counts") %>%
  # na.omit() %>%
  identity()

pseudobulk_DE_23 = run_de(expr_input,
                          meta = seu_23@meta.data,
                          de_family = "pseudobulk",
                          de_method = "DESeq2",
                          de_type = "LRT",
                          min_reps = 2,
                          min_cells = 2)


inspect_pseudobulk_DE_23 <-
  pseudobulk_DE_23 %>%
  inspect_libra_output()




