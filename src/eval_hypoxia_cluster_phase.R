#!/usr/bin/env Rscript
# Evaluate hypoxia-flagged clusters in cell-cycle-phase terms.
#
# Motivation: the cluster-based hypoxia split (split_hypoxia_by_clusters) moves
# whole clusters into the high-hypoxia subset purely on a mean-hypoxia_score
# outlier rule, independent of hypoxia marker-gene presence. That can pull in
# clusters that are really *cell-cycle phase* states (e.g. histone/PLK2/EZH2 or
# heat-shock clusters) rather than an orthogonal hypoxia axis. Filtering such a
# cluster would delete a fine-grained phase state.
#
# For every cluster that is ever flagged into the high subset (across
# resolutions 0.2/0.4/0.6, single round) this script quantifies:
#   (1) phase purity  - dominant-phase fraction + Shannon entropy over G1/S/G2M
#   (2) CC score profile - mean S.Score / G2M.Score of the flagged cluster
#   (3) cross-resolution phase divergence - whether a res-0.2 flagged cluster
#       fragments into res-0.6 subclusters with different dominant phases
#   (4) hypoxia-vs-cycle attribution - paired with n_hyp_markers / top_markers
#
# The hypsplit_r1_res.* cluster columns are NOT persisted in the labeled .rds,
# so the split's clustering is re-derived here with identical parameters
# (deterministic; n_iter = 1 clusters all cells in round 1).
#
# Outputs (results/hypoxia_cluster_split/):
#   hypoxia_flagged_cluster_phase_eval.csv   one row per flagged (sample,res,cluster)
#   hypoxia_flagged_crossres_detail.csv      res0.2 cluster -> res0.6 subcluster
#   hypoxia_flagged_crossres_summary.csv     phase-splitting summary per res0.2 cluster
#   hypoxia_flagged_cluster_phase_eval.pdf   overview plots (best-effort)

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat)
  library(dplyr); library(readr); library(stringr); library(glue)
})

RESOLUTIONS <- c(0.2, 0.6, 1.0)
PHASES      <- c("G1", "S", "G2M")
PURE_FRAC   <- 0.70   # dominant-phase fraction >= this => "phase-restricted"
HYP_MIN     <- 2L     # n_hyp_markers >= this => hypoxia-marker-positive

seu_dir <- "output/seurat"
out_dir <- "results/hypoxia_cluster_split"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional single-/subset-sample run: HYP_EVAL_SAMPLES=SRX...,SRX... restricts the
# sample set; HYP_EVAL_SUFFIX appends to output filenames so a quick per-sample
# run cannot clobber the full run's CSVs.
SUFFIX      <- Sys.getenv("HYP_EVAL_SUFFIX", "")
SAMPLE_ENV  <- Sys.getenv("HYP_EVAL_SAMPLES", "")
outp <- function(base, ext) file.path(out_dir, paste0(base, SUFFIX, ext))

shannon_bits <- function(p) { p <- p[p > 0]; if (length(p) == 0) 0 else -sum(p * log2(p)) }

# Samples: every sample that produced a split log (i.e. was ever clustered).
log_all_path <- file.path(out_dir, "hypoxia_split_log_all.csv")
if (file.exists(log_all_path)) {
  samples <- sort(unique(readr::read_csv(log_all_path, show_col_types = FALSE)$sample_id))
} else {
  samples <- stringr::str_extract(
    list.files(seu_dir, pattern = "_seu_hypoxia_labeled\\.rds$"), "SR[RX][0-9]+")
  samples <- sort(unique(samples[!is.na(samples)]))
}
if (nzchar(SAMPLE_ENV)) {
  want <- trimws(strsplit(SAMPLE_ENV, ",")[[1]])
  samples <- intersect(samples, want)
}
message("Evaluating ", length(samples), " samples",
        if (nzchar(SUFFIX)) paste0(" (suffix '", SUFFIX, "')") else "", ".")

rows <- list()   # per flagged cluster
xres <- list()   # cross-resolution detail

for (sid in samples) {
  lp <- file.path(seu_dir, paste0(sid, "_seu_hypoxia_labeled.rds"))
  if (!file.exists(lp)) { message("skip ", sid, ": no labeled rds"); next }
  message("== ", sid, " ==")
  seu <- readRDS(lp)

  da <- if ("gene" %in% names(seu@assays)) "gene" else DefaultAssay(seu)
  DefaultAssay(seu) <- da
  if (!"data" %in% SeuratObject::Layers(seu[[da]]))
    seu <- Seurat::NormalizeData(seu, assay = da, verbose = FALSE)

  # Ensure cell-cycle Phase / scores (reuse the project's 1q-free scoring).
  if (!all(c("Phase", "S.Score", "G2M.Score") %in% colnames(seu@meta.data))) {
    seu <- tryCatch(annotate_cell_cycle_without_1q(seu, organism = "human"),
                    error = function(e) {
                      message("  CC scoring failed: ", conditionMessage(e)); seu })
  }
  has_phase <- "Phase" %in% colnames(seu@meta.data)
  has_hyp   <- "hypoxia_score" %in% colnames(seu@meta.data)

  if (!"pca" %in% names(seu@reductions)) { message("  no pca; skip"); rm(seu); gc(); next }

  # Re-derive the split's clustering (identical params to split_hypoxia_by_clusters).
  n_cells <- ncol(seu); k_param <- min(20L, n_cells - 1L)
  seu <- Seurat::FindNeighbors(seu, dims = 1:30, reduction = "pca",
                               graph.name = paste0(da, c("_nn", "_snn")),
                               k.param = k_param, verbose = FALSE)

  for (res in RESOLUTIONS) {
    grp <- paste0("hypsplit_r1_res.", res)
    seu <- Seurat::FindClusters(seu, graph.name = paste0(da, "_snn"),
                                resolution = res, verbose = FALSE)
    seu@meta.data[[grp]] <- seu$seurat_clusters
    det <- identify_hypoxia_clusters(seu, grp, mad_k = 3, return_detail = TRUE)
    if (!is.data.frame(det) || nrow(det) == 0) next
    flagged <- as.character(det$cluster[det$flagged])
    for (cl in flagged) {
      inx <- as.character(seu@meta.data[[grp]]) == cl
      n   <- sum(inx, na.rm = TRUE); if (n == 0) next
      d   <- det[as.character(det$cluster) == cl, , drop = FALSE][1, ]
      if (has_phase) {
        tab  <- table(factor(seu$Phase[inx], levels = PHASES))
        frac <- as.numeric(tab) / sum(tab)
        dom  <- PHASES[which.max(tab)]; domf <- max(frac); ent <- shannon_bits(frac)
      } else { frac <- rep(NA_real_, 3); dom <- NA; domf <- NA; ent <- NA }
      rows[[length(rows) + 1]] <- data.frame(
        sample_id = sid, resolution = res, cluster = cl, n_cells = n,
        frac_G1 = frac[1], frac_S = frac[2], frac_G2M = frac[3],
        dominant_phase = dom, dominant_frac = domf, phase_entropy_bits = ent,
        mean_S   = if (has_phase) mean(seu$S.Score[inx],   na.rm = TRUE) else NA_real_,
        mean_G2M = if (has_phase) mean(seu$G2M.Score[inx], na.rm = TRUE) else NA_real_,
        mean_hypoxia_score = if (has_hyp) mean(seu$hypoxia_score[inx], na.rm = TRUE) else NA_real_,
        n_hyp_markers = d$n_hyp_markers, matched_genes = d$matched_genes,
        top_markers = d$top_markers, cluster_mean_score = d$mean_score,
        outlier_fence = d$outlier_fence,
        stringsAsFactors = FALSE)
    }
  }

  # Cross-resolution: how each flagged cluster at the LOWEST resolution
  # partitions across the HIGHEST resolution's clusters.
  res_lo <- min(RESOLUTIONS); res_hi <- max(RESOLUTIONS)
  c_lo <- paste0("hypsplit_r1_res.", res_lo); c_hi <- paste0("hypsplit_r1_res.", res_hi)
  if (all(c(c_lo, c_hi) %in% colnames(seu@meta.data)) && has_phase) {
    det_lo <- identify_hypoxia_clusters(seu, c_lo, mad_k = 3, return_detail = TRUE)
    det_hi <- identify_hypoxia_clusters(seu, c_hi, mad_k = 3, return_detail = TRUE)
    fl_lo  <- as.character(det_lo$cluster[det_lo$flagged])
    fl_hi  <- as.character(det_hi$cluster[det_hi$flagged])
    for (cl in fl_lo) {
      inx <- as.character(seu@meta.data[[c_lo]]) == cl
      th  <- sort(table(as.character(seu@meta.data[[c_hi]][inx])), decreasing = TRUE)
      for (sc in names(th)) {
        cells <- inx & as.character(seu@meta.data[[c_hi]]) == sc
        pht   <- table(factor(seu$Phase[cells], levels = PHASES))
        xres[[length(xres) + 1]] <- data.frame(
          sample_id = sid, res_lo = res_lo, res_hi = res_hi,
          res_lo_cluster = cl, res_hi_subcluster = sc,
          n_cells = as.integer(th[[sc]]),
          frac_of_parent = as.numeric(th[[sc]]) / sum(th),
          res_hi_flagged = sc %in% fl_hi,
          dominant_phase = PHASES[which.max(pht)],
          dominant_frac  = if (sum(pht) > 0) max(pht) / sum(pht) else NA_real_,
          stringsAsFactors = FALSE)
      }
    }
  }
  rm(seu); gc()
}

per_cluster <- dplyr::bind_rows(rows)
if (nrow(per_cluster) > 0) {
  per_cluster <- per_cluster |>
    dplyr::mutate(
      phase_restricted = !is.na(dominant_frac) & dominant_frac >= PURE_FRAC,
      hypoxia_marker_positive = n_hyp_markers >= HYP_MIN,
      # Interpretation: safe-to-filter == genuine hypoxia, phase-mixed.
      call = dplyr::case_when(
        !hypoxia_marker_positive &  phase_restricted ~ "phase_artifact",
        !hypoxia_marker_positive & !phase_restricted ~ "ambiguous_nonhypoxia",
         hypoxia_marker_positive &  phase_restricted ~ "hypoxia_but_phase_restricted",
        TRUE                                          ~ "genuine_hypoxia_phase_mixed"))
  readr::write_csv(per_cluster, outp("hypoxia_flagged_cluster_phase_eval", ".csv"))
  message("wrote per-cluster eval: ", nrow(per_cluster), " flagged clusters")
} else message("no flagged clusters found")

xr <- dplyr::bind_rows(xres)
if (nrow(xr) > 0) {
  readr::write_csv(xr, outp("hypoxia_flagged_crossres_detail", ".csv"))
  xr_sum <- xr |>
    dplyr::filter(frac_of_parent >= 0.10) |>          # ignore tiny fragments
    dplyr::group_by(sample_id, res_lo, res_hi, res_lo_cluster) |>
    dplyr::summarise(
      n_sub_hi = dplyr::n_distinct(res_hi_subcluster),
      n_distinct_dominant_phase = dplyr::n_distinct(dominant_phase),
      dominant_phases = paste(sort(unique(dominant_phase)), collapse = "/"),
      phase_splitting = dplyr::n_distinct(dominant_phase) > 1,
      .groups = "drop")
  readr::write_csv(xr_sum, outp("hypoxia_flagged_crossres_summary", ".csv"))
  message("wrote cross-resolution summary: ", nrow(xr_sum), " res-0.2 flagged clusters")
}

# ---- overview plots (best-effort) ----
tryCatch({
  suppressPackageStartupMessages({ library(ggplot2); library(tidyr); library(patchwork) })
  if (exists("per_cluster") && nrow(per_cluster) > 0) {
    pc <- per_cluster |> dplyr::mutate(id = paste(sample_id, resolution, cluster, sep = ":"))
    long <- pc |>
      dplyr::select(id, sample_id, call, frac_G1, frac_S, frac_G2M) |>
      tidyr::pivot_longer(c(frac_G1, frac_S, frac_G2M),
                          names_to = "phase", values_to = "frac") |>
      dplyr::mutate(phase = factor(sub("frac_", "", phase), levels = PHASES))
    p1 <- ggplot(long, aes(id, frac, fill = phase)) +
      geom_col() + coord_flip() +
      facet_grid(sample_id ~ ., scales = "free_y", space = "free_y") +
      labs(title = "Phase composition of hypoxia-flagged clusters",
           x = NULL, y = "fraction of cells", fill = "Phase") +
      theme_bw(base_size = 7) +
      theme(strip.text.y = element_text(angle = 0), legend.position = "top")
    p2 <- ggplot(pc, aes(mean_S, mean_G2M, colour = call, size = n_cells)) +
      geom_point(alpha = 0.8) +
      labs(title = "Flagged clusters: cell-cycle vs hypoxia attribution",
           x = "mean S.Score", y = "mean G2M.Score") +
      theme_bw(base_size = 9) + theme(legend.position = "right")
    ggsave(outp("hypoxia_flagged_cluster_phase_eval", ".pdf"),
           p1 / p2 + patchwork::plot_layout(heights = c(3, 1)),
           width = 10, height = max(8, 0.18 * nrow(pc) + 4), limitsize = FALSE)
    message("wrote overview pdf")
  }
}, error = function(e) message("plotting failed: ", conditionMessage(e)))

cat("DONE\n")
