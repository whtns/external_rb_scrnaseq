#!/usr/bin/env Rscript
# The panel-agreement question, redone on COORDINATES instead of labels.
#
# Earlier passes joined heatmap and bulk segment sets on the `seg` label and
# reported 20/39 samples disagreeing. That was an artefact: numbat re-derives
# segment labels per stage (155 renamed between segs_consensus and joint_post,
# 58 between it and bulk_clones), and by coordinate all 688 called segments are
# present in BOTH fields, with zero truly absent.
#
# So the open question is only whether the two panels' FILTERS diverge --
# heatmap needs a cell with p_cnv >= p_min, bulk needs a clone with
# LLR >= min_LLR -- once the segments are matched properly.
#
# Read-only.

suppressPackageStartupMessages({ library(dplyr) })
P_MIN <- 0.9; MIN_LLR <- 5
fld <- function(nb, n) tryCatch(nb[[n]], error = function(e) NULL)

best_over <- function(chrom, start, end, tbl, value_col) {
  cand <- tbl[tbl$CHROM == chrom, , drop = FALSE]
  if (!nrow(cand)) return(NA_real_)
  ov <- pmin(end, cand$seg_end) - pmax(start, cand$seg_start)
  ln <- pmin(end - start, cand$seg_end - cand$seg_start)
  keep <- ov > 0 & (ov / pmax(ln, 1)) >= 0.9
  if (!any(keep, na.rm = TRUE)) return(NA_real_)
  suppressWarnings(max(cand[[value_col]][keep], na.rm = TRUE))
}

rows <- lapply(sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds")), function(f) {
  s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
  nb <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(nb)) return(NULL)
  segs <- fld(nb,"segs_consensus"); jp <- fld(nb,"joint_post"); bc <- fld(nb,"bulk_clones")
  if (is.null(segs)||is.null(jp)||is.null(bc)) return(NULL)
  need <- c("CHROM","seg_start","seg_end")
  if (!all(need %in% names(segs)) || !all(need %in% names(jp)) || !all(need %in% names(bc))) return(NULL)
  sc <- if ("cnv_state_post" %in% names(segs)) "cnv_state_post" else "cnv_state"
  sg <- segs |> filter(.data[[sc]] != "neu") |> distinct(seg, .keep_all = TRUE)
  if (!nrow(sg)) return(NULL)
  purrr::map_dfr(seq_len(nrow(sg)), function(i) {
    r <- sg[i,]
    tibble::tibble(sample_id = s, seg = r$seg,
      max_p_cnv = best_over(r$CHROM, r$seg_start, r$seg_end, jp, "p_cnv"),
      max_LLR   = best_over(r$CHROM, r$seg_start, r$seg_end, bc, "LLR"))
  })
})
d <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
d <- d |> mutate(in_hm = !is.na(max_p_cnv) & max_p_cnv >= P_MIN,
                 in_bk = !is.na(max_LLR)   & max_LLR   >= MIN_LLR,
                 mismatch = in_hm != in_bk)
per <- d |> group_by(sample_id) |> summarise(n = n(), mm = sum(mismatch), .groups="drop")
cat(sprintf("=== coordinate-matched, at the CURRENT defaults (p_min %.2f, min_LLR %g) ===\n", P_MIN, MIN_LLR))
cat(sprintf("segments: %d   mismatched: %d (%.1f%%)\n", nrow(d), sum(d$mismatch), 100*mean(d$mismatch)))
cat(sprintf("samples agreeing exactly: %d / %d\n", sum(per$mm == 0), nrow(per)))
if (any(per$mm > 0)) {
  cat("\nsamples still disagreeing:\n"); print(as.data.frame(per |> filter(mm>0)), row.names=FALSE)
  cat("\nthe disagreeing segments:\n")
  print(as.data.frame(d |> filter(mismatch) |> select(sample_id, seg, max_p_cnv, max_LLR, in_hm, in_bk)), row.names=FALSE)
}
readr::write_csv(d, "results/panel_agreement_by_coord.csv")
cat("\nwrote results/panel_agreement_by_coord.csv\nDIAG DONE\n")
