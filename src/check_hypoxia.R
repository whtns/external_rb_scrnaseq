#!/usr/bin/env Rscript

source("packages.R")
source("functions.R")

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)

tar_load("cluster_orders")

# seus <- c(
#     "SRR13884249" = "output/seurat/SRR13884249_filtered_seu.rds",
# 	"SRR14800534" = "output/seurat/SRR14800534_filtered_seu.rds",
# 	"SRR14800535" = "output/seurat/SRR14800535_filtered_seu.rds",
#     "SRR14800536" = "output/seurat/SRR14800536_filtered_seu.rds",
#     "SRR13884246" = "output/seurat/SRR13884246_filtered_seu.rds",
#     "SRR17960484" = "output/seurat/SRR17960484_filtered_seu.rds"
# ) |> 
# 	set_names()
# 
# seus <- set_names(seus, str_extract(names(seus), "SRR[0-9]*"))
# 
# seus <- seus |> 
# 	map(readRDS)
# 
# 
# subset_seu_by_clones <- function(seu, clones){
# 	seu <- seu[,seu$clone_opt %in% clones]
# 	
# 	return(seu)
# }
# 
# # integrated 1q_16q seu ------------------------------
# 
# seus <- map2(seus, list(
# 	"SRR13884249" = c(1,2),
# 	"SRR14800534" = c(1,2,3),
# 	"SRR14800535" = c(1,2,3),
# 	"SRR14800536" = c(1,2,3),
# 	"SRR13884246" = c(1,2),
# 	"SRR17960484" = c(1,2)
# 	),
# 	subset_seu_by_clones)
# 

# integrated_seu <- seuratTools::integration_workflow(seus)
# 
# integrated_seu <- ScaleData(integrated_seu)
# 
# integrated_seu <- seurat_reduce_dimensions(integrated_seu)
# 
# integrated_seu <- seurat_cluster(integrated_seu, resolution = seq(0.2, 2.0, by = 0.2), seurat_assay = "integrated")
# 
# integrated_seu <- add_hypoxia_score(integrated_seu)
# 
# 
# saveRDS(integrated_seu, "output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6.rds")

integrated_seu <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6.rds")

plot_hypoxia_score(integrated_seu) + 
    geom_hline(yintercept = 0.4, color = "red") + 
    geom_hline(yintercept = 0.5, color = "red")

# ## hypoxia < 0.5 ------------------------------
# 
# hypoxia_less_0_5_seu <- integrated_seu[,integrated_seu$hypoxia_score < 0.5] |> 
#     seurat_cluster(resolution = seq(0.4, 0.8, by = 0.2), seurat_assay = "integrated")
# 
rds_path <- "output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6_less_0_5.rds"

# saveRDS(hypoxia_less_0_5_seu, rds_path)

hypoxia_less_0_5_seu <- readRDS(rds_path)
# 
# ## hypoxia >= 0.5 ------------------------------
# 
# hypoxia_more_0_5_seu <- integrated_seu[,integrated_seu$hypoxia_score >= 0.5] |> 
#     seurat_cluster(resolution = seq(0.4, 0.8, by = 0.2), seurat_assay = "integrated")
# 
rds_path <- "output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6_more_0_5.rds"

# saveRDS(hypoxia_more_0_5_seu, rds_path)

hypoxia_more_0_5_seu <- readRDS(rds_path)
# 
# ## hypoxia < 0.4 ------------------------------
# 
# hypoxia_less_0_4_seu <- integrated_seu[,integrated_seu$hypoxia_score < 0.4] |> 
#     seurat_cluster(resolution = seq(0.4, 0.8, by = 0.2), seurat_assay = "integrated")
# 
rds_path <- "output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6_less_0_4.rds"

# saveRDS(hypoxia_less_0_4_seu, rds_path)

hypoxia_less_0_4_seu <- readRDS(rds_path)
# 
# ## hypoxia >= 0.4 ------------------------------
# 
# hypoxia_more_0_4_seu <- integrated_seu[,integrated_seu$hypoxia_score >= 0.4] |> 
#     seurat_cluster(resolution = seq(0.4, 0.8, by = 0.2), seurat_assay = "integrated")
# 
rds_path <- "output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6_more_0_4.rds"

# saveRDS(hypoxia_more_0_4_seu, rds_path)

hypoxia_more_0_4_seu <- readRDS(rds_path)

# 1q from six 1q/16q ------------------------------

integrated_seu <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6.rds")

# debug(select_1q_clones)
# debug(subset_to_1q)
# debug(make_clone_distribution_figure_debug)

pdf1 <- select_1q_clones(integrated_seu,
                         file_id = "integrated_seu_1q_16q_afterall6.rds",
                         slug = "_1q",
                         # group.by = "integrated_snn_res.0.8",
                         cluster_orders = cluster_orders)

pdf2 <- select_1q_clones(hypoxia_less_0_5_seu,
                         file_id = "integrated_seu_1q_16q_afterall6_less_0_5.rds",
                         slug = "_1q_less_0_5",
                         # group.by = "integrated_snn_res.0.8",
                         cluster_orders = cluster_orders)

pdf3 <- select_1q_clones(hypoxia_more_0_5_seu,
                         file_id = "integrated_seu_1q_16q_afterall6_more_0_5.rds",
                         slug = "_1q_more_0_5",
                         # group.by = "integrated_snn_res.0.8",
                         cluster_orders = cluster_orders)

pdf4 <- select_1q_clones(hypoxia_less_0_4_seu,
                         file_id = "integrated_seu_1q_16q_afterall6_less_0_4.rds",
                         slug = "_1q_less_0_4",
                         # group.by = "integrated_snn_res.0.8",
                         cluster_orders = cluster_orders)

pdf5 <- select_1q_clones(hypoxia_more_0_4_seu,
                         file_id = "integrated_seu_1q_16q_afterall6_more_0_4.rds",
                         slug = "_1q_more_0_4",
                         # group.by = "integrated_snn_res.0.8",
                         cluster_orders = cluster_orders)

# 16q for only three samples from six 1q/16q ------------------------------

integrated_seu <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall6.rds")

# debug(select_16q_clones)
# # debug(subset_to_1q)
# debug(make_clone_distribution_figure_debug)
# debug(plot_distribution_of_clones_across_clusters)

pdf6 <- select_16q_clones(integrated_seu, 
                          file_id = "integrated_seu_1q_16q_afterall6.rds", 
                          slug = "_16q", 
                          # group.by = "integrated_snn_res.0.8", 
                          cluster_orders = cluster_orders)

pdf7 <- select_16q_clones(hypoxia_less_0_5_seu, 
                          file_id = "integrated_seu_1q_16q_afterall6_less_0_5.rds", 
                          slug = "_16q_less_0_5", 
                          # group.by = "integrated_snn_res.0.8",
                          cluster_orders = cluster_orders)

pdf8 <- select_16q_clones(hypoxia_more_0_5_seu, 
                          file_id = "integrated_seu_1q_16q_afterall6_more_0_5.rds", 
                          slug = "_16q_more_0_5", 
                          # group.by = "integrated_snn_res.0.8",
                          cluster_orders = cluster_orders)

pdf9 <- select_16q_clones(hypoxia_less_0_4_seu, 
                          file_id = "integrated_seu_1q_16q_afterall6_less_0_4.rds", 
                          slug = "_16q_less_0_4", 
                          # group.by = "integrated_snn_res.0.8",
                          cluster_orders = cluster_orders)

pdf10 <- select_16q_clones(hypoxia_more_0_4_seu, 
                           file_id = "integrated_seu_1q_16q_afterall6_more_0_4.rds", 
                           slug = "_16q_more_0_4", 
                           # group.by = "integrated_snn_res.0.8",
                           cluster_orders = cluster_orders)

qpdf::pdf_combine(c(pdf1, pdf2, pdf3, pdf4, pdf5, pdf6, pdf7, pdf8, pdf9, pdf10), "results/hypoxia_test2.pdf")


# # trio seu ------------------------------
# 
# seus <- seus[c("SRR14800534", "SRR14800535", "SRR14800536")]
# 
# trio_seus <- map2(seus, list(
#     "SRR14800534" = c(1,2,3),
#     "SRR14800535" = c(1,2,3),
#     "SRR14800536" = c(1,2,3)
# ),
# subset_seu_by_clones)
# 
# integrated_seu <- seuratTools::integration_workflow(trio_seus)
# 
# integrated_seu <- ScaleData(integrated_seu)
# 
# integrated_seu <- seurat_reduce_dimensions(integrated_seu)
# 
# integrated_seu <- seurat_cluster(integrated_seu, resolution = seq(0.2, 2.0, by = 0.2), seurat_assay = "integrated")
# 
# saveRDS(integrated_seu, "output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall_trio.rds")
# 
# integrated_seu <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_16q_afterall_trio.rds")
# 
# 
# # integrated 16q seu ------------------------------
# 
# clone_selection <- tribble(
#     ~batch, ~clone_opt,
#     "SRR14800534", c(1,2),
#     "SRR14800535", c(1,2),
#     "SRR14800536", c(1,2)
# ) |> 
#     tidyr::unnest(clone_opt)
# 
# selected_cells <- integrated_seu@meta.data |> 
#     rownames_to_column("cell") |> 
#     dplyr::inner_join(clone_selection, by = c("batch", "clone_opt"))
# 
# seu_16q <- integrated_seu[,selected_cells$cell]
# 
# seu_16q$scna <-
#     factor(ifelse(str_detect(seu_16q$GT_opt, "16"), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))
# 
# saveRDS(seu_16q, "output/seurat/integrated_1q_16q/integrated_seu_16q_afterall.rds")
# 
# seu_16q |> 
#     SplitObject(split.by = "batch") |> 
#     imap(~saveRDS(.x, glue("output/seurat/integrated_1q_16q/{.y}_integrated_16q_filtered_seu.rds")))
