#!/usr/bin/env Rscript
# Render plot_seu_marker_heatmap for SRX10264523_hypoxia_low with the auto phase
# labels wired in; copy the resulting PDF to a stable path for inspection.
suppressPackageStartupMessages({
  source("packages.R")   # same package set the pipeline attaches
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")
})
dir.create("tmp", showWarnings = FALSE)

nb_paths <- list.files("output/numbat_sridhar", pattern = "_numbat\\.rds$", full.names = TRUE)
clone_simpl <- yaml::read_yaml("config/large_clone_simplifications.yaml")

p <- plot_seu_marker_heatmap(
  "output/seurat/SRX10264523_hypoxia_low_seu.rds",
  cluster_order = NULL,        # -> dummy_cluster_order -> resolution SCT_snn_res.0.6
  nb_paths = nb_paths,
  clone_simplifications = clone_simpl,
  tmp_plot_path = TRUE
)
cat("returned:", p, "\n")
if (!is.na(p) && file.exists(p)) {
  out <- "results/TEST_SRX10264523_autophase.pdf"
  file.copy(p, out, overwrite = TRUE)
  cat("copied to:", out, "\n")
}
cat("DONE\n")
