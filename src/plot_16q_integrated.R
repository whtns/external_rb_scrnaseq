source("packages.R")
source("functions.R")

library(speckle)


seu_path <- "output/seurat/SRX11133594_SRX11133593_SRX11133592_seu.rds"

seu <- readRDS(seu_path)

seu <- seu[,!is.na(seu$clone_opt)]

seu <- seu[,seu$clone_opt %in% c(1,2,3)]

DefaultAssay(seu) <- "gene"

sce <- as.SingleCellExperiment(seu)

# debug(propeller)
# debug(propeller.anova)

propeller(sce, clusters = sce$abbreviation, sample=sce$batch, group=sce$clone_opt, transform = "asin")


test0 <- make_integrated_numbat_plots(seu_path)

sample_id = "SRX11133594_SRX11133593_SRX11133592"

test0$markerplot
ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_sample_marker.pdf"), height = 4, width = 5)

test0$dimplot
ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_dimplot.pdf"), width = 8, height = 4)

test0$distplot
ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_clone_distribution.pdf"), width = 4, height = 4)
