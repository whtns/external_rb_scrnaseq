#!/usr/bin/env Rscript
# Find a (p_min, min_LLR) pair that makes the heatmap and the bulk-clone panel
# draw exactly the same segments (github #42 follow-up, #11).
#
# The two panels filter the SAME non-neutral segs_consensus set, but on
# different statistics:
#
#   heatmap     keeps a segment if ANY CELL has  p_cnv >= p_min      (default 0.9)
#   bulk panel  keeps a segment if ANY CLONE has LLR   >= min_LLR    (default 5)
#
# Those are different quantities -- a per-cell posterior and a pseudobulk
# log-likelihood ratio -- so there is no a priori reason a threshold pair exists
# that equalises them. This measures whether one does.
#
# Pass 1 reads each RDS once and caches, per (sample, segment), the max p_cnv
# over cells and the max LLR over clones. Pass 2 sweeps thresholds over that
# cached table, which is cheap.
#
# Read-only with respect to the numbat objects and the targets store.

suppressPackageStartupMessages({ library(dplyr) })

CACHE <- "results/panel_segment_stats.csv"

fld <- function(nb, n) tryCatch(nb[[n]], error = function(e) NULL)

# ---- pass 1: per-segment statistics -----------------------------------------
if (!file.exists(CACHE)) {
  rds <- sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds"))
  cat("pass 1: reading", length(rds), "objects\n")

  rows <- lapply(rds, function(f) {
    s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
    nb <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(nb)) return(NULL)
    segs <- fld(nb, "segs_consensus"); jp <- fld(nb, "joint_post"); bc <- fld(nb, "bulk_clones")
    if (is.null(segs) || is.null(jp) || is.null(bc)) return(NULL)

    state_col <- if ("cnv_state_post" %in% names(segs)) "cnv_state_post" else "cnv_state"
    called <- segs$seg[segs[[state_col]] != "neu"]
    if (!length(called)) return(NULL)

    p <- jp |> filter(seg %in% called) |> group_by(seg) |>
      summarise(max_p_cnv = max(p_cnv, na.rm = TRUE), .groups = "drop")
    l <- bc |> filter(seg %in% called, !is.na(LLR)) |> group_by(seg) |>
      summarise(max_LLR = max(LLR, na.rm = TRUE), .groups = "drop")

    tibble::tibble(sample_id = s, seg = called) |>
      left_join(p, by = "seg") |> left_join(l, by = "seg")
  })

  d <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
  readr::write_csv(d, CACHE)
  cat("wrote", CACHE, "-", nrow(d), "sample-segments\n\n")
} else {
  cat("pass 1: reusing", CACHE, "\n\n")
}

d <- readr::read_csv(CACHE, show_col_types = FALSE) |>
  mutate(max_p_cnv = ifelse(is.na(max_p_cnv), -Inf, max_p_cnv),
         max_LLR   = ifelse(is.na(max_LLR),   -Inf, max_LLR))

n_samples <- dplyr::n_distinct(d$sample_id)
cat("segments:", nrow(d), " samples:", n_samples, "\n")

# ---- pass 2: sweep -----------------------------------------------------------
score <- function(p_min, min_llr) {
  x <- d |>
    mutate(in_hm = max_p_cnv >= p_min, in_bulk = max_LLR >= min_llr,
           mismatch = in_hm != in_bulk)
  per <- x |> group_by(sample_id) |>
    summarise(mm = sum(mismatch), n_drawn = sum(in_hm | in_bulk), .groups = "drop")
  tibble::tibble(
    p_min = p_min, min_LLR = min_llr,
    samples_agreeing = sum(per$mm == 0),
    segs_mismatched  = sum(per$mm),
    segs_drawn       = sum(per$n_drawn))
}

p_grid <- c(0, 0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
l_grid <- c(0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50)

grid <- do.call(rbind, lapply(p_grid, function(p)
  do.call(rbind, lapply(l_grid, function(l) score(p, l)))))

cat("\n=== current defaults ===\n")
print(as.data.frame(grid |> filter(p_min == 0.9, min_LLR == 5)), row.names = FALSE)

exact <- grid |> filter(samples_agreeing == n_samples)

cat("\n=== threshold pairs giving EXACT agreement on all", n_samples, "samples ===\n")
if (nrow(exact) == 0) {
  cat("  none in the swept grid\n")
} else {
  # among exact pairs, prefer the one that still DRAWS the most segments -- an
  # exact match achieved by drawing nothing is not useful.
  print(as.data.frame(exact |> arrange(desc(segs_drawn), desc(p_min), desc(min_LLR))),
        row.names = FALSE)
  best <- exact |> arrange(desc(segs_drawn), desc(p_min), desc(min_LLR)) |> slice(1)
  cat(sprintf("\nrecommended: p_min = %g, min_LLR = %g  (draws %d segments)\n",
              best$p_min, best$min_LLR, best$segs_drawn))
}

cat("\n=== best 12 pairs by agreement, then by segments drawn ===\n")
print(as.data.frame(grid |> arrange(desc(samples_agreeing), desc(segs_drawn)) |> head(12)),
      row.names = FALSE)

readr::write_csv(grid, "results/panel_threshold_sweep.csv")
cat("\nwrote results/panel_threshold_sweep.csv\nSWEEP DONE\n")
