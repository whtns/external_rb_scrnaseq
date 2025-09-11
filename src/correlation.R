#!/usr/bin/env Rscript   

library(targets)
suppressPackageStartupMessages(source("packages.R"))
source("functions.R")
library(plotgardener)
library(tidymodels)

seu <- readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall.rds")

# correlation ------------------------------

seurat.data = seu[["integrated"]]@data
seurat.data.cor.pearson = cor(t(as.matrix(seurat.data)), 
                              method = "pearson")
# write.csv(seurat.data.cor, file =
# '../results/2019.12.16/E16.5_hipp_Seurat_correlations.csv'
# )

hypoxia_genes <-
    dplyr::filter(seu@misc$markers$clusters$presto, Cluster == "hypoxia_2") |> 
    slice_head(n = 10) |>
    dplyr::pull(Gene.Name) |>
    identity()
    
mt_genes <- str_subset(rownames(seu), "MT-.*")
# mt_genes <- c("MT-CO3", "MT-ND3", "MT-CYB", "MT-ATP6", "MT-CO2")

hypoxia_genes <- c(hypoxia_genes, mt_genes)


partial.coex.pearson = seurat.data.cor.pearson[rownames(seurat.data.cor.pearson) %in% 
                                                   hypoxia_genes, colnames(seurat.data.cor.pearson) %in% 
                                                   hypoxia_genes]
diag(partial.coex.pearson) = 0

seurat.data.cor.spearman = cor(t(as.matrix(seurat.data)), 
                               method = "spearman")
# write.csv(seurat.data.cor, file =
# '../results/2019.12.16/E16.5_hipp_Seurat_correlations.csv'
# )

partial.coex.spearman = seurat.data.cor.spearman[rownames(seurat.data.cor.spearman) %in% 
                                                     hypoxia_genes, colnames(seurat.data.cor.spearman) %in% 
                                                     hypoxia_genes]
diag(partial.coex.spearman) = 0

# pearson correlation ------------------------------

partial.coex.pearson = reshape2::melt(partial.coex.pearson)
colnames(partial.coex.pearson) = c("g1","g2","corr")

partial.coex.pearson$g1 <- factor(partial.coex.pearson$g1, hypoxia_genes)
partial.coex.pearson$g2 <- factor(partial.coex.pearson$g2, hypoxia_genes)

P = ggplot(partial.coex.pearson) + 
    geom_tile(aes(x=g1,y=g2, fill = corr),colour = "black", show.legend = TRUE) +
    #  facet_grid( g1 ~ g2  ,scales = "free", space = "free") + 
    scale_fill_gradient2(mid = "white",limits=c(-1, 1),low = "#DC0000B2", high = "#3C5488B2")+
    #scale_fill_gradient2(low = "darkred", mid = "white",  high = "darkblue", midpoint = 0,na.value = "grey80", space = "Lab", guide = "colourbar", aesthetics = "fill", limits = lim_coex, oob=scales::squish)+ theme(legend.position="bottom")+
    theme(#legend.title = element_blank(),
        #strip.text.x = element_text(color = "red"),
        #axis.text.y = element_text(color = ),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1.0),
        legend.position="bottom"
    ) #+geom_text(aes(label=ifelse(t_hk == "hk", "H","")), color="grey", size=3)
P

ggsave("results/pearson_corr.pdf", w = 6, h = 6) |> 
    browseURL()

# spearman ------------------------------

partial.coex.spearman = reshape2::melt(partial.coex.spearman)
colnames(partial.coex.spearman) = c("g1","g2","corr")

partial.coex.spearman$g1 <- factor(partial.coex.spearman$g1, hypoxia_genes)
partial.coex.spearman$g2 <- factor(partial.coex.spearman$g2, hypoxia_genes)

S = ggplot(partial.coex.spearman) + 
    geom_tile(aes(x=g1,y=g2, fill = corr),colour = "black", show.legend = TRUE) +
    #  facet_grid( g1 ~ g2  ,scales = "free", space = "free") + 
    scale_fill_gradient2(mid = "white",limits=c(-1, 1),low = "#DC0000B2", high = "#3C5488B2")+
    #scale_fill_gradient2(low = "darkred", mid = "white",  high = "darkblue", midpoint = 0,na.value = "grey80", space = "Lab", guide = "colourbar", aesthetics = "fill", limits = lim_coex, oob=scales::squish)+ theme(legend.position="bottom")+
    theme(#legend.title = element_blank(),
        #strip.text.x = element_text(color = "red"),
        #axis.text.y = element_text(color = ),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1.0),
        legend.position="bottom"
    ) #+geom_text(aes(label=ifelse(t_hk == "hk", "H","")), color="grey", size=3)
S

ggsave("results/spearman_corr.pdf", w = 6, h = 6) |> 
    browseURL()
