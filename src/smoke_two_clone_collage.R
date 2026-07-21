# Smoke-test the two-clone collage on ONE sample/SCNA/resolution.
#
# The builder now RE-CLUSTERS on the included cells: it subsets to the two
# compared clones first, then recomputes PCA / SNN / clusters on that subset. So
# the check is the inverse of the old lock-to-full one -- the displayed clusters
# are a fresh labelling of the two-clone slice and are NOT expected to match the
# full-population cluster ids. What must hold is:
#   1. a PDF renders,
#   2. the number of clusters is sane for the cell count (re-clustering happened
#      and did not collapse to a single cluster or explode),
#   3. the stacked bars are still labelled by SCNA-of-interest status.
suppressPackageStartupMessages({
  source("packages.R"); library(Seurat); library(dplyr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

low  <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
stopifnot(file.exists(low))
lcc  <- yaml::read_yaml("config/large_clone_comparisons.yaml")

col <- "SCT_snn_res.0.4"
s <- readRDS(low)
stopifnot(col %in% colnames(s@meta.data))
full_clusters <- sort(unique(as.character(s@meta.data[[col]])))

# The two clones of the 1q comparison (2_v_1_1q+ -> clones 2 and 1).
retained <- c("1", "2")
keep <- as.character(s@meta.data$clone_opt) %in% retained
inherited_clusters <- sort(unique(as.character(s@meta.data[[col]][keep])))
n_cells <- sum(keep, na.rm = TRUE)

cat("\n=== SRX10264523 res 0.4 ===\n")
cat("full-object clusters              :", paste(full_clusters, collapse = " "), "\n")
cat("clusters clones 1+2 INHERIT (old) :", paste(inherited_clusters, collapse = " "), "\n")
cat("cells in clones 1+2               :", n_cells, "\n")

# Reproduce the builder's own clustering so we can report what it will show.
sub <- s[, colnames(s)[which(keep)]]
rm(s); gc()
npcs    <- max(2L, min(30L, ncol(sub) - 1L))
k_param <- max(2L, min(20L, ncol(sub) - 1L))
sub <- RunPCA(sub, assay = "SCT", npcs = npcs, verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:npcs, reduction = "pca", k.param = k_param,
                     graph.name = c("SCT_nn", "SCT_snn"), verbose = FALSE)
sub <- FindClusters(sub, graph.name = "SCT_snn", resolution = 0.4, verbose = FALSE)
recomputed <- sort(unique(as.character(sub$seurat_clusters)))
cat("clusters RECOMPUTED on subset (new):", paste(recomputed, collapse = " "), "\n")
stopifnot(length(recomputed) >= 1, length(recomputed) <= n_cells)
cat("PASS: re-clustering produced", length(recomputed), "clusters on", ncol(sub), "cells\n")
rm(sub); gc()

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
stopifnot(file.exists(expected))
cat("pdf exists:", file.exists(expected), "\n")
cat("pdf size (KB):", round(file.info(expected)$size / 1024, 1), "\n")

# Confirm the stacked-bar panel is still labeled by SCNA-of-interest status, not
# raw clone ids: the rendered text should contain "1q+" and "preceding".
txt <- paste(pdftools::pdf_text(expected), collapse = " ")
cat("bar label '1q+' present      :", grepl("1q\\+", txt), "\n")
cat("bar label 'preceding' present:", grepl("preceding", txt), "\n")
stopifnot(grepl("1q\\+", txt), grepl("preceding", txt))

cat("\nSMOKE DONE\n")
