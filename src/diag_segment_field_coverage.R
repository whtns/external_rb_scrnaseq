#!/usr/bin/env Rscript
# Why do joint_post and bulk_clones cover different segment sets? (github #42/#11)
#
# 99 of 688 unique called segments exist in one field of the numbat object but
# not the other -- 98 in bulk_clones only, 1 in joint_post only. That asymmetry
# is upstream of any threshold, so tuning p_min / min_LLR to make the panels
# agree would paper over it rather than fix it.
#
# This asks what distinguishes a segment that reaches joint_post from one that
# does not. Candidate explanations, all testable from the object:
#
#   * it is an evidence threshold  -> the absent ones have low LLR / few
#     genes / few SNPs
#   * it is the CNV state          -> particular states never reach joint_post
#   * it is size                   -> tiny segments are dropped
#   * it is a naming mismatch      -> the segment is there under another label,
#     i.e. the seg sets are differently indexed rather than genuinely different
#
# Read-only.

suppressPackageStartupMessages({ library(dplyr) })

fld <- function(nb, n) tryCatch(nb[[n]], error = function(e) NULL)

rds <- sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) rds <- rds[grepl(paste(args, collapse = "|"), rds)]

rows <- lapply(rds, function(f) {
  s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) return(NULL)
  segs <- fld(nb, "segs_consensus"); jp <- fld(nb, "joint_post"); bc <- fld(nb, "bulk_clones")
  if (is.null(segs) || is.null(jp) || is.null(bc)) return(NULL)

  state_col <- if ("cnv_state_post" %in% names(segs)) "cnv_state_post" else "cnv_state"
  sg <- segs |> filter(.data[[state_col]] != "neu") |> distinct(seg, .keep_all = TRUE)
  if (!nrow(sg)) return(NULL)

  in_jp <- sg$seg %in% unique(jp$seg)
  in_bc <- sg$seg %in% unique(bc$seg)

  # per-segment attributes that might explain the split
  grab <- function(col) if (col %in% names(sg)) sg[[col]] else rep(NA_real_, nrow(sg))

  tibble::tibble(
    sample_id = s,
    seg       = sg$seg,
    state     = sg[[state_col]],
    in_joint_post  = in_jp,
    in_bulk_clones = in_bc,
    seg_length = grab("seg_length"),
    n_genes    = grab("n_genes"),
    n_snps     = grab("n_snps"),
    LLR        = grab("LLR"),
    LLR_x      = grab("LLR_x"),
    LLR_y      = grab("LLR_y")
  )
})

d <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
stopifnot(nrow(d) > 0)

d <- d |> mutate(group = case_when(
  in_joint_post &  in_bulk_clones ~ "both",
  !in_joint_post &  in_bulk_clones ~ "bulk_clones only",
  in_joint_post & !in_bulk_clones ~ "joint_post only",
  TRUE ~ "neither"))

cat("=== how many called segments reach each field ===\n")
print(as.data.frame(d |> count(group, name = "n_segments")), row.names = FALSE)

cat("\n=== do they differ by CNV state? ===\n")
print(as.data.frame(d |> count(group, state, name = "n") |>
        tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)),
      row.names = FALSE)

cat("\n=== do they differ by evidence / size? (median per group) ===\n")
print(as.data.frame(d |> group_by(group) |> summarise(
  n          = n(),
  seg_length = median(seg_length, na.rm = TRUE),
  n_genes    = median(n_genes, na.rm = TRUE),
  n_snps     = median(n_snps, na.rm = TRUE),
  LLR        = median(LLR, na.rm = TRUE),
  LLR_x      = median(LLR_x, na.rm = TRUE),
  LLR_y      = median(LLR_y, na.rm = TRUE),
  .groups = "drop")), row.names = FALSE)

# Is it simply that numbat keeps only segments clearing its own min_LLR (5) in
# the per-cell step? If so, "bulk_clones only" should sit almost entirely below
# that line and "both" almost entirely above.
cat("\n=== LLR vs numbat's own min_LLR = 5 ===\n")
print(as.data.frame(d |> filter(group %in% c("both", "bulk_clones only")) |>
  group_by(group) |> summarise(
    n = n(),
    below_5 = sum(LLR < 5, na.rm = TRUE),
    pct_below_5 = round(100 * mean(LLR < 5, na.rm = TRUE), 1),
    .groups = "drop")), row.names = FALSE)

cat("\n=== the joint_post-only segments (should be rare) ===\n")
jo <- d |> filter(group == "joint_post only")
if (nrow(jo)) print(as.data.frame(jo |> select(sample_id, seg, state, LLR, n_genes, n_snps)),
                    row.names = FALSE) else cat("  none\n")

readr::write_csv(d, "results/segment_field_coverage.csv")
cat("\nwrote results/segment_field_coverage.csv\nDIAG DONE\n")
