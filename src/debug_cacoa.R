#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(cacoa)

seu <- readRDS("output/seurat/CopyOfSRR14800534_SRR14800535_SRR14800536_seu.rds")

seu <- seu[,!is.na(seu$clone_opt)]



# sample.groups: vector with condition labels per sample named with sample ids
# cell.groups: cell type annotation vector named by cell ids
# sample.per.cell: vector with sample labels per cell named with cell ids
# ref.level: id of the condition, corresponding to the reference (i.e. control)
# target.level: id of the condition, corresponding to the target (i.e. case)

cao <- cacoa::Cacoa$new(
  seu,
  sample.groups=c(),
  sample.per.cell=seu$batch,
  cell.groups=seu$abbreviation,
  target.level='IPF',
  ref.level='Control',
  n.cores=4,
  verbose=FALSE
)

# set default plot parameters
cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))
cao$plot.theme <- cao$plot.theme + theme(legend.background=element_blank())
