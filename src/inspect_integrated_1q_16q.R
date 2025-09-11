#!/usr/bin/env Rscript

source("packages.R")
source("functions.R")

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)

seu_16q <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_16q_afterall.rds")

seu_1q <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall.rds")

# dimplots ------------------------------
DimPlot(seu_16q, group.by = "integrated_snn_res.0.4") + 
    DimPlot(seu_16q, group.by = "Phase", split.by = "scna")

DimPlot(seu_1q, group.by = "integrated_snn_res.0.4") + 
    DimPlot(seu_1q, group.by = "Phase", split.by = "scna")

# bar plot ------------------------------

seu_meta <- seu_16q@meta.data %>%
    identity()

var_x = "scna"
var_y = "integrated_snn_res.0.4"

summarized_clones <-
    seu_meta %>%
    dplyr::select(.data[[var_x]], .data[[var_y]]) %>%
    dplyr::mutate(scna = "all")

ggplot(seu_meta) +
    geom_bar(position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]]), color = "black", width = 0.6) +
    geom_bar(data = summarized_clones, position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]]), color = "black", width = 0.6) +
    # scale_x_discrete(limits = rev) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20), limits = rev) +
    coord_flip()
