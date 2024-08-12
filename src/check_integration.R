library(seuratTools)
library(tidyverse)

seu <- seu[,seu$batch %in% c("SRR13884247", "SRR13884248", "SRR17960484")]

seu$batch <- factor(seu$batch)

seus <- Seurat::SplitObject(seu, split.by = "batch")

seus |> 
	integration_workflow(find_markers = FALSE) |> 
	saveRDS("output/seurat/integrated_2p/seurat_2p_integrated_trio.rds")

source("packages.R")
source("functions.R")

# debug(make_clone_distribution_figure)

seu <- find_all_markers(seu, seurat_assay = "integrated")

seu <- find_all_markers(seu, metavar = "integrated_snn_res.0.2")

seu <- find_all_markers(seu, metavar = "integrated_snn_res.0.4")

seu <- find_all_markers(seu, metavar = "integrated_snn_res.0.6")

saveRDS(seu, "output/seurat/integrated_2p/seurat_2p_integrated_trio.rds")

test0 <- make_clone_distribution_figure("output/seurat/integrated_2p/seurat_2p_integrated_trio.rds", scna_of_interest = "2p\\+", width =20, height = 10, group.bys = c("integrated_snn_res.0.2", "integrated_snn_res.0.4", "integrated_snn_res.0.6"))

cc_data <-FetchData(seu, c("G2M.Score", "S.Score", paste0("integrated_snn_res.0.", c(2,4,6)), "scna", "batch"))

ggplot(cc_data, aes(x = S.Score, y = G2M.Score, color = batch)) +
	geom_point() + 
	facet_wrap(~`integrated_snn_res.0.6`)

seus <- Seurat::SplitObject(seu, split.by = "batch")

imap(seus, ~saveRDS(.x, glue("output/seurat/integrated_2p/{.y}_filtered_seu.rds")))

?seurat_find_markers()

seu_list
seu0 <- find_all_markers(seu, metavar = "clusters")

?find_all_markers