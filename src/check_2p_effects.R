source("packages.R")
source("functions.R")
library(targets)

library(tidyverse)

tar_load('debranched_seus_2p')

# cluster_terms <- map(debranched_seus_2p, enrich_by_cluster)


tar_load(c("numbat_rds_files", "large_clone_comparisons"))

debug(make_clone_comparison)

diffex_scna_enriched_clusters <- function(seu_path, scna_clusters = c("g1_0", "g1_4"), phase = "g1", scna_of_interest = "2p") {
  
  seu <- readRDS(seu_path)
  
  g1_seu <- seu[,seu$phase_level == phase]
  
  g1_seu$divergence = ifelse(g1_seu$clusters %in% scna_clusters, scna_of_interest, "diploid")
  
  test0 <- FindMarkers(g1_seu, ident.1 = scna_of_interest, ident.2 = "diploid", group.by = "divergence")
  
  test1 <- enrichment_analysis(test0, gene_set = "hallmark")
  
  plot_enrichment(test1)

}

diffex_scna_enriched_clusters(debranched_seus_2p[[1]], scna_clusters = c("g1_0", "g1_4"))

diffex_scna_enriched_clusters(debranched_seus_2p[[2]], scna_clusters = c("g1_2", "g1_8"))

diffex_scna_enriched_clusters(debranched_seus_2p[[6]], scna_clusters = c("g1_0"))

g1_2p_effects <- map(debranched_seus_2p[c(1,2,4,5,6)], find_diffex_clones_in_phase, phase = "g1", scna_of_interest = "2p", numbat_rds_files, large_clone_comparisons, location = "all")

pdf("~/results/tmp.pdf")
print(test0)
dev.off()

browseURL("~/results/tmp.pdf")
