# Smoke-test the lock-to-full two-clone collage on ONE sample/SCNA/resolution.
# Confirms: (1) plot_scna_two_clone_res_collages renders a PDF for SRX10264523 1q
# at res 0.4, and (2) the clusters it displays are a SUBSET of the full-population
# clusters at the SAME persisted column (i.e. a true zoom, not a re-clustering).
suppressPackageStartupMessages({
  source("packages.R"); library(Seurat); library(dplyr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

low  <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
stopifnot(file.exists(low))
lcc  <- yaml::read_yaml("config/large_clone_comparisons.yaml")

# Full-object cluster set at the persisted column we will zoom into.
s <- readRDS(low)
col <- "SCT_snn_res.0.4"
stopifnot(col %in% colnames(s@meta.data))
full_clusters <- sort(unique(as.character(s@meta.data[[col]])))

# The two clones of the 1q comparison (2_v_1_1q+ -> clones 2 and 1).
retained <- c("1", "2")
keep <- as.character(s@meta.data$clone_opt) %in% retained
two_clone_clusters <- sort(unique(as.character(s@meta.data[[col]][keep])))
cat("\n=== SRX10264523 SCT_snn_res.0.4 ===\n")
cat("full-object clusters      :", paste(full_clusters, collapse = " "), "\n")
cat("clones 1+2 occupy clusters:", paste(two_clone_clusters, collapse = " "), "\n")
cat("subset is subset of full  :", all(two_clone_clusters %in% full_clusters), "\n")
cat("cells in clones 1+2       :", sum(keep, na.rm = TRUE), "\n")
rm(s); gc()

# Render just res 0.4 via the real function.
out <- plot_scna_two_clone_res_collages(
  seu_path                = low,
  scna_of_interest        = "1q",
  large_clone_comparisons = lcc,
  resolutions             = 0.4,
  nb_paths                = NULL
)
cat("\nreturned:", paste(out, collapse = ", "), "\n")
expected <- "results/SRX10264523_hypoxia_low_seu.rds__1q_res0.4_heatmap_phase_scatter_patchwork.pdf"
cat("pdf exists:", file.exists(expected), "\n")
if (file.exists(expected))
  cat("pdf size (KB):", round(file.info(expected)$size / 1024, 1), "\n")
cat("\nSMOKE DONE\n")
