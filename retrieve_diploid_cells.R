setwd("/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj")
suppressPackageStartupMessages(source("packages.R"))
lapply(list.files("R", full.names = TRUE, pattern = "[.]R$"), source) |>
  invisible()

diploid_seu <- targets::tar_read(diploid_seu, store = "_targets_r431")
diploid_seu
