#!/usr/bin/env Rscript
# Render plot_seu_marker_heatmap for SRX11133594_hypoxia_low with clones (clone_opt)
# displayed in place of scna; copy the resulting PDF to a stable path for inspection.
suppressPackageStartupMessages({
  source("packages.R")   # same package set the pipeline attaches
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")
})
dir.create("tmp", showWarnings = FALSE)

nb_paths <- list.files("output/numbat_sridhar", pattern = "_numbat\\.rds$", full.names = TRUE)
clone_simpl <- yaml::read_yaml("config/large_clone_simplifications.yaml")

p <- plot_seu_marker_heatmap(
  "output/seurat/SRX11133594_hypoxia_low_seu.rds",
  cluster_order = NULL,
  nb_paths = nb_paths,
  clone_simplifications = clone_simpl,
  tmp_plot_path = TRUE
)
cat("returned:", p, "\n")
if (!is.na(p) && file.exists(p)) {
  out <- "results/TEST_SRX11133594_clone_swap.pdf"
  file.copy(p, out, overwrite = TRUE)
  cat("copied to:", out, "\n")
}
cat("DONE\n")
