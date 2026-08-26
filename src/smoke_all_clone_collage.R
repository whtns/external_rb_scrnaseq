# Smoke test for plot_scna_all_clone_res_collages(): one 1q sample, two
# resolutions, exercising the hypoxia-removed clustree marking + removed-cell
# UMAP band before committing to a full rebuild.
suppressMessages({
  # Mirror the pipeline's attached-package environment exactly (find_all_markers
  # from seuratTools, igraph's as.igraph, ggraph's edge_colourbar guide, clustree,
  # ...) rather than hand-picking libraries -- packages.R is what the successful
  # two-clone rebuild had attached. Then overlay the latest numbatHelpers source.
  source("packages.R")
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")
})
tar_config_set(store = "_targets_r431")

sid  <- "SRX10831287"
seu  <- file.path("output/seurat", paste0(sid, "_filtered_seu.rds"))
high <- file.path("output/seurat", paste0(sid, "_hypoxia_high_seu.rds"))
stopifnot(file.exists(seu), file.exists(high))

nb_paths <- tar_read(numbat_rds_files)
clone_simp <- tryCatch(tar_read(large_clone_simplifications), error = function(e) NULL)

cat("=== running plot_scna_all_clone_res_collages on", sid, "===\n")
out <- plot_scna_all_clone_res_collages(
  seu,
  scna_of_interest   = "1q",
  resolutions        = c(0.4, 0.8),
  high_hypoxia_paths = high,
  nb_paths           = nb_paths,
  clone_simplifications = clone_simp
)

cat("\n=== returned paths ===\n"); print(out)
ok <- !is.na(out) & file.exists(out)
cat("PDFs written:", sum(ok), "/", length(out), "\n")
for (p in out[ok]) cat("  ", p, "  (", round(file.size(p)/1024), "KB )\n")
if (sum(ok) == 0) quit(status = 1)
cat("SMOKE PASS\n")
