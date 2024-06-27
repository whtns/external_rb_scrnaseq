#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(patchwork)
library(Seurat)
library(ggpubr)

seu <- readRDS("output/seurat/SRR14800534_filtered_seu.rds")

# total ------------------------------
total_genes <- c("MT1X", "MT2A", "CENPF", "CKS1B")

seu_total <-
  seu %>%
  FetchData(c(total_genes, "scna")) %>%
  dplyr::mutate(scna = factor(scna))


# cluster ------------------------------
cluster_genes <- c("KIF14", "COPA", "EPRS")

seu_cluster <-
  seu %>%
  subset(subset = SCT_snn_res.0.6 == 6) %>%
  FetchData(c(cluster_genes, "scna")) %>%
  dplyr::mutate(scna = factor(scna))

mycomparisons = list(c("", "16q-"), c("16q-", "16q- 1q+"))

myvplot <- function(mygene, seu){
  ggpubr::ggviolin(seu, x = "scna", y = mygene, fill = "scna",
                   add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(label.y = y_lim) +
    stat_compare_means(comparisons = mycomparisons, label = "p.signif") +
    scale_y_continuous(limits = c(-1, y_lim+1), breaks = seq(0, y_lim+1, step)) +
    guides(fill = "none") +
    NULL

}

# per cluster ------------------------------

y_lim = 3
step = 2

map(cluster_genes, myvplot, seu_cluster) %>%
  wrap_plots(ncol = 1) +
  plot_annotation(title = "1q+ differentially expressed")

ggsave("~/tmp/Rplot.pdf", height = 6)

browseURL("~/tmp/Rplot.pdf")

# total ------------------------------

y_lim = 7
step = 2

map(total_genes, myvplot, seu_total) %>%
  wrap_plots(ncol = 1) +
  # plot_annotation(title = "1q+ differentially expressed") +
  NULL

ggsave("~/tmp/Rplot.pdf", height = 8)

browseURL("~/tmp/Rplot.pdf")
