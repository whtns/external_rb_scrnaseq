#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(SingleCellExperiment)
library(Seurat)
library(glue)

rename_assay_to_gene <- function(seu){
  seu@assays$gene <- seu@assays$originalexp

  DefaultAssay(seu) <- "gene"

  return(seu)
}

processed_seus <- dir_ls("data/ALSF/", glob = "*processed.rds", recurse = TRUE) %>%
  set_names(str_extract(., "SCPCL[0-9]*")) %>%
  map(readRDS) %>%
  map(as.Seurat, counts = "counts", data = "logcounts") %>%
  map(rename_assay_to_gene) %>%
  identity()

imap(processed_seus, ~saveRDS(.x, glue("output/seurat/{.y}_processed_seu.rds")))

merged_seu <- purrr::reduce(processed_seus, merge)

get_ensgene_from_symbol <- function(mysymbol){
  annotables::grch38 %>%
    dplyr::filter(symbol %in% mysymbol) %>%
    dplyr::pull(ensgene)
}

test0 <- map(processed_seus, ~FeaturePlot(.x, features = get_ensgene_from_symbol("TFF1")))

symbols_to_rownames <- function(sce){

  sce <- processed_seus$SCPCL000752

  rd <- rowData(sce)

  rd <- rd[!is.na(rd$gene_symbol),]

  sce <- sce[rownames(rd),]

  symbol_counts <- counts(sce)
  rownames(symbol_counts) <- rd$gene_symbol

  seu <- CreateSeuratObject(counts = symbol_counts, meta.data = as.data.frame(colData(sce)))

  # Set feature metadata, AKA rowData. Super intuitive, right?
  rownames(rd) <- rd$gene_symbol
  seu[["RNA"]]@meta.features <- as.data.frame(rd)

  return(minimal_seu)

}

