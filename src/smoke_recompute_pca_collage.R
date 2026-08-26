# Smoke test the recompute_pca + capped-sweep (0.2..0.8) low-hypoxia collage path.
# Uses SRX10264526 (lost the most cells at split -> best test of the fresh-PCA
# re-clustering) and reports, per resolution, the produced PDF size AND the
# resulting cluster-size distribution so we can confirm fragmentation dropped.
suppressPackageStartupMessages({
  source("packages.R")
  library(targets); library(Seurat); library(dplyr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")

p <- "output/seurat/SRX10264526_hypoxia_low_seu.rds"
res <- seq(0.2, 0.8, by = 0.2)

out <- plot_hypoxia_low_res_collages(
  p, nb_paths = NULL, resolutions = res, recompute_pca = TRUE)

cat("\n===== produced collages =====\n")
for (f in out) {
  if (is.na(f)) { cat("  -> NA\n"); next }
  sz <- if (file.exists(f)) file.info(f)$size else NA
  cat(sprintf("  %s | %s bytes | %s\n", basename(f), format(sz, big.mark=","),
              if (!is.na(sz) && sz > 20000) "OK" else "!! TOO SMALL"))
}

# Independently recompute PCA + clusters the same way and report cluster sizes,
# to confirm the small-cluster burden vs the persisted (inherited-PCA) columns.
cat("\n===== cluster sizes: fresh-PCA re-clustering (SRX10264526) =====\n")
s <- readRDS(p); DefaultAssay(s) <- "SCT"
s <- RunPCA(s, assay = "SCT", npcs = 30, verbose = FALSE)
s <- FindNeighbors(s, dims = 1:30, reduction = "pca",
                   graph.name = c("SCT_nn","SCT_snn"), verbose = FALSE)
for (r in res) {
  s <- FindClusters(s, graph.name = "SCT_snn", resolution = r, verbose = FALSE)
  tab <- sort(table(s$seurat_clusters), decreasing = TRUE)
  cat(sprintf("  res %.1f: %d clusters | min %d | median %.0f | n<20 = %d\n",
              r, length(tab), min(tab), median(tab), sum(tab < 20)))
}
cat("\nSMOKE DONE\n")
