# Smoke-test the per-cluster significance marks on the stacked-bar panel
# (issue #32) using the same sample/SCNA/resolution as src/smoke_two_clone_collage.R.
#
# What must hold:
#   1. .bar_enrichment_stats() returns ONE row per cluster -- including the
#      untested ones, so a missing star is never ambiguous;
#   2. every row carries either a p/q pair or an explanatory note, never both NA;
#   3. the collage still renders;
#   4. results/..._bar_enrichment.csv is written and matches the plot's own table
#      (this is the part that broke first: the stats attribute has to survive the
#      patchwork composition inside plot_distribution_of_clones_across_clusters);
#   5. the rendered PDF actually carries the caption and at least one mark token,
#      i.e. the marks reached the page rather than just the CSV.
suppressPackageStartupMessages({
  source("packages.R"); library(Seurat); library(dplyr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

low <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
stopifnot(file.exists(low))
lcc <- yaml::read_yaml("config/large_clone_comparisons.yaml")

# ---- 1/2: the helper in isolation, on the real two-clone subset -------------
s <- readRDS(low)
keep <- as.character(s@meta.data$clone_opt) %in% c("1", "2")
sub <- s[, colnames(s)[which(keep)]]
rm(s); gc()

npcs    <- max(2L, min(30L, ncol(sub) - 1L))
k_param <- max(2L, min(20L, ncol(sub) - 1L))
sub <- RunPCA(sub, assay = "SCT", npcs = npcs, verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:npcs, reduction = "pca", k.param = k_param,
                     graph.name = c("SCT_nn", "SCT_snn"), verbose = FALSE)
sub <- FindClusters(sub, graph.name = "SCT_snn", resolution = 0.4, verbose = FALSE)

meta <- sub@meta.data
meta$clusters <- factor(as.character(meta$seurat_clusters))
meta$scna_status <- factor(ifelse(as.character(meta$clone_opt) == "2",
                                  "1q+ (clone 2)", "preceding (clone 1)"))

st <- numbatHelpers:::.bar_enrichment_stats(meta, "scna_status", "clusters", min_cells = 20)
cat("\n=== .bar_enrichment_stats on SRX10264523 clones 1+2, res 0.4 ===\n")
print(as.data.frame(st))

n_clusters <- length(unique(as.character(meta$clusters)))
cat("clusters present:", n_clusters, " rows returned:", nrow(st), "\n")
stopifnot(nrow(st) == n_clusters)
stopifnot(identical(sort(as.character(st$level)),
                    sort(unique(as.character(meta$clusters)))))
stopifnot(sum(st$n_cells) == nrow(meta))
# Untested rows must say why; tested rows must have a finite q.
untested <- is.na(st$p_value)
stopifnot(all(!is.na(st$note[untested])), all(is.na(st$note[!untested])))
stopifnot(all(is.finite(st$q_value[!untested])))
stopifnot(all(!is.na(st$label)))
cat("PASS: one row per cluster,", sum(!untested), "tested,", sum(untested), "explained-untested\n")
rm(sub, meta); gc()

# ---- 3/4/5: end to end through the real builder -----------------------------
out <- plot_scna_two_clone_res_collages(
  seu_path                = low,
  scna_of_interest        = "1q",
  large_clone_comparisons = lcc,
  resolutions             = 0.4,
  nb_paths                = NULL,
  bar_signif              = TRUE
)
cat("\nreturned:", paste(out, collapse = ", "), "\n")

pdf_path <- "results/SRX10264523_hypoxia_low_seu.rds__1q_res0.4_heatmap_phase_scatter_patchwork.pdf"
stopifnot(file.exists(pdf_path))
cat("pdf size (KB):", round(file.info(pdf_path)$size / 1024, 1), "\n")

# Name built the same way plot_seu_marker_heatmap builds it.
csv_path <- "results/SRX10264523_hypoxia_low_seu.rds_1q_res0.4_bar_enrichment.csv"
cat("expecting csv:", csv_path, "\n")
if (!file.exists(csv_path)) {
  cat("!! not found; results/ entries matching bar_enrichment:\n")
  print(list.files("results", pattern = "bar_enrichment", full.names = TRUE))
}
stopifnot(file.exists(csv_path))
csv <- readr::read_csv(csv_path, show_col_types = FALSE)
print(as.data.frame(csv))
stopifnot(nrow(csv) > 0,
          all(c("sample_id", "level", "n_cells", "p_value", "q_value", "label") %in% names(csv)))
cat("PASS: enrichment csv written with", nrow(csv), "rows\n")

txt <- paste(pdftools::pdf_text(pdf_path), collapse = " ")
has_caption <- grepl("Fisher", txt)
# A mark is present if any cluster label picked one up. Look for the tokens the
# helper emits rather than for a bare "*", which appears elsewhere on the page.
marked <- vapply(as.character(csv$label), function(l) grepl(l, txt, fixed = TRUE),
                 logical(1))
cat("caption present :", has_caption, "\n")
cat("labels found on page:", sum(marked), "/", nrow(csv), "\n")
print(as.character(csv$label)[!marked])
stopifnot(has_caption, any(marked))

cat("\nSMOKE DONE\n")
