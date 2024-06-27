library(chevreul)
library(tidyverse)
library(glue)
library(sceasy)
library(chevreul)
# devtools::load_all()
matplotlib <<- reticulate::import("matplotlib", convert = TRUE)
matplotlib$use("Agg", force = TRUE)
pyplot <<- reticulate::import("matplotlib.pyplot", delay_load = TRUE)
scvelo <<- reticulate::import("scvelo", delay_load = TRUE)

seu_path <- "~/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat/SRR14800534_filtered_seu.rds"

seu <- readRDS(seu_path)

# seu <- seu[rownames(seu) %in% VariableFeatures(seu$SCT),]

seu@reductions$tsne@cell.embeddings[,1] <- seu$S.Score
seu@reductions$tsne@cell.embeddings[,2] <- seu$G2M.Score

sample_id <- str_extract(seu_path, "SRR[0-9]*")

debug(prep_scvelo)
debug(run_scvelo)
file.remove(glue("/home/skevin/single_cell_projects/resources/external_rb_scrnaseq_proj/output/velocyto/{sample_id}.h5ad"))

test0 <- prep_scvelo(seu, glue("/home/skevin/single_cell_projects/resources/external_rb_scrnaseq_proj/output/velocyto/{sample_id}.loom"), velocity_mode = "deterministic", assay = "gene")

scvelo$pl$velocity_embedding_grid(test0, save = glue("{sample_id}_grid.png"), color = "SCT_snn_res.0.6", dpi = 200, figsize = c(20, 12), basis = "tsne")

scvelo$pl$velocity_embedding_grid(test0, save = glue("{sample_id}_grid_phase.png"), color = "Phase", dpi = 200, figsize = c(20, 12), basis = "tsne",title = sample_id)

scvelo$pl$velocity_embedding(test0, basis = "umap", color = "SCT_snn_res.0.6", arrow_length = 4, arrow_size = 1, dpi = 200, figsize = c(20, 12), save= glue("{sample_id}.png"))

plot_scvelo(test0, group.by = "SCT_snn_res.0.6", plot_method = "arrow", save= glue("{sample_id}.png"), basis = "tsne")

plot_scvelo(test0, group.by = "Phase", plot_method = "arrow", save= glue("{sample_id}_phase.png"), basis = "tsne")

scvelo$pl$velocity_embedding_grid(test0, save = glue("{sample_id}_grid_phase.png"), color = "Phase", dpi = 200, figsize = c(20, 12), basis = "tsne", title = sample_id)


# scvelo$pl$velocity_embedding_grid(test0, save = glue("{sample_id}_grid_scna.png"), color = "scna", dpi = 200, figsize = c(20, 12), basis = "tsne")
#
# scvelo$pl$velocity(test0, c('ARL6IP1',  'RRM2', 'PTTG1', 'TFF1'), ncols=2, save = "asdf.png", basis="tsne")
