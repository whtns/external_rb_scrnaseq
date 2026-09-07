#!/usr/bin/env Rscript
# Do the sample-summary panels actually show the SAME SEGMENTS? (github #42, #11)
#
# Fixing #42 guaranteed provenance: the heatmap and the bulk-clone panel now
# read the same *_numbat.rds at the same selected round. That is NOT the same as
# guaranteeing they display the same segments, because they draw different
# fields of that object, over different cell sets, through different filters:
#
#   heatmap      nb$joint_post   per-cell posterior, Seurat cell subset, p_min 0.9
#   bulk panel   nb$bulk_clones  pseudobulk per clone, all numbat cells, min_LLR 5
#   segs         nb$segs_consensus   the segmentation both are derived from
#
# This reports, per sample, the segment sets each panel can draw and where they
# disagree. Read-only.

suppressPackageStartupMessages({ library(dplyr) })

MIN_LLR <- 5      # plot_psbulk default
P_MIN   <- 0.9    # make_numbat_heatmaps default, as used by the summary targets

rds <- sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) rds <- rds[grepl(paste(args, collapse = "|"), rds)]

fld <- function(nb, n) tryCatch(nb[[n]], error = function(e) NULL)

rows <- lapply(rds, function(f) {
  s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) return(NULL)

  segs <- fld(nb, "segs_consensus")
  jp   <- fld(nb, "joint_post")
  bc   <- fld(nb, "bulk_clones")
  if (is.null(segs) || is.null(jp) || is.null(bc)) return(NULL)

  state_col <- if ("cnv_state_post" %in% names(segs)) "cnv_state_post" else "cnv_state"
  seg_all <- segs$seg[segs[[state_col]] != "neu"]

  # what the heatmap can draw: non-neutral segments with at least one cell over
  # the posterior threshold
  jp_hit <- jp |>
    filter(seg %in% seg_all, p_cnv >= P_MIN) |>
    distinct(seg) |> pull(seg)

  # what the bulk panel can draw: segments clearing the pseudobulk LLR in at
  # least one clone
  bc_hit <- bc |>
    filter(seg %in% seg_all, !is.na(LLR), LLR >= MIN_LLR) |>
    distinct(seg) |> pull(seg)

  data.frame(
    sample_id      = s,
    n_segs_called  = length(seg_all),
    n_heatmap      = length(jp_hit),
    n_bulk         = length(bc_hit),
    n_shared       = length(intersect(jp_hit, bc_hit)),
    heatmap_only   = paste(sort(setdiff(jp_hit, bc_hit)), collapse = ";"),
    bulk_only      = paste(sort(setdiff(bc_hit, jp_hit)), collapse = ";"),
    stringsAsFactors = FALSE)
})

d <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
stopifnot(nrow(d) > 0)

d$agree <- d$n_heatmap == d$n_shared & d$n_bulk == d$n_shared

cat("=== segments each panel can draw, per sample ===\n")
print(as.data.frame(d[, c("sample_id","n_segs_called","n_heatmap","n_bulk","n_shared","agree")]),
      row.names = FALSE)

cat(sprintf("\nsamples where the two panels draw the SAME segment set: %d / %d\n",
            sum(d$agree), nrow(d)))

if (any(!d$agree)) {
  cat("\n=== disagreements ===\n")
  for (i in which(!d$agree)) {
    cat(sprintf("%s  (called %d)\n", d$sample_id[i], d$n_segs_called[i]))
    if (nzchar(d$heatmap_only[i])) cat("    heatmap only: ", d$heatmap_only[i], "\n", sep = "")
    if (nzchar(d$bulk_only[i]))    cat("    bulk only   : ", d$bulk_only[i], "\n", sep = "")
  }
}

readr::write_csv(d, "results/panel_segment_agreement.csv")
cat("\nwrote results/panel_segment_agreement.csv\nDIAG DONE\n")
