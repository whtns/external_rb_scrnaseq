#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

# excel_path <- "http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"
#
# download.file(excel_path, "data/table_S3_surfaceome.xlsx")

sheets <- readxl::excel_sheets("data/table_S3_surfaceome.xlsx") %>% set_names(.)
surface_markers <- map(sheets, ~readxl::read_excel("data/table_S3_surfaceome.xlsx", .x, skip = 1)) %>%
  map(janitor::clean_names)

# wtc11 expression
download.file("https://s3-us-west-2.amazonaws.com/downloads.allencell.org/genome-sequence/AICStranscriptomics_v1.2.zip", "data/wtc11_expression_allen.zip")

wtc11_tpm <- read_csv("data/wtc11_expression_allen/AICStranscriptomics032018/march2018_tpm.csv")

markers_16q <-
  surface_markers$SurfaceomeMasterTable %>%
  dplyr::left_join(wtc11_tpm, by = c("uni_prot_gene" = "gene")) %>%
  dplyr::left_join(annotables::grch38, by = c("uni_prot_gene" = "symbol"), relationship = "many-to-many") %>%
  dplyr:::filter(chr == "16") %>%
  dplyr::filter(surfaceome_label_source == "pos. trainingset")
