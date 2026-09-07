#!/usr/bin/env Rscript
# Are the "missing" segments actually missing, or just differently LABELLED?
# (github #42 / #11 follow-up)
#
# The previous diagnostic ruled out every evidence-based explanation: segments
# absent from joint_post have HIGHER median LLR (9.7 vs 9.3) than those present,
# none fall below numbat's own min_LLR of 5, and they are larger with more SNPs.
# The 57 absent from BOTH fields have a median LLR of 62.7 -- the strongest in
# the cohort. Strong segments do not get silently filtered; strong segments get
# silently RENAMED.
#
# numbat re-derives segment labels ("1a", "1b", "16d", ...) per stage, and
# segs_consensus has been through relevel_chrom(). If joint_post and
# bulk_clones carry their own labelling, then matching the panels on `seg`
# compares labels, not genomic intervals.
#
# This tests that directly: for each segment absent from a field by NAME, is
# there a segment in that field covering the same CHROM/start/end?
#
# Read-only.

suppressPackageStartupMessages({ library(dplyr) })

`%||%` <- function(a, b) if (is.null(a)) b else a

fld <- function(nb, n) tryCatch(nb[[n]], error = function(e) NULL)

rds <- sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) rds <- rds[grepl(paste(args, collapse = "|"), rds)]

# coordinate-level match: same chromosome, and intervals overlapping by >=90%
# of the shorter one (labels shift, boundaries move slightly between stages)
overlaps <- function(chrom, start, end, tbl) {
  if (is.null(tbl) || !all(c("CHROM","seg_start","seg_end") %in% names(tbl))) return(NA)
  cand <- tbl[tbl$CHROM == chrom, , drop = FALSE]
  if (!nrow(cand)) return(FALSE)
  ov <- pmin(end, cand$seg_end) - pmax(start, cand$seg_start)
  ln <- pmin(end - start, cand$seg_end - cand$seg_start)
  any(ov > 0 & (ov / pmax(ln, 1)) >= 0.9, na.rm = TRUE)
}

rows <- lapply(rds, function(f) {
  s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) return(NULL)
  segs <- fld(nb, "segs_consensus"); jp <- fld(nb, "joint_post"); bc <- fld(nb, "bulk_clones")
  if (is.null(segs) || is.null(jp) || is.null(bc)) return(NULL)
  if (!all(c("CHROM","seg_start","seg_end") %in% names(segs))) return(NULL)

  state_col <- if ("cnv_state_post" %in% names(segs)) "cnv_state_post" else "cnv_state"
  sg <- segs |> filter(.data[[state_col]] != "neu") |> distinct(seg, .keep_all = TRUE)
  if (!nrow(sg)) return(NULL)

  jp_u <- if (all(c("CHROM","seg_start","seg_end") %in% names(jp)))
            distinct(jp, CHROM, seg, seg_start, seg_end) else NULL
  bc_u <- if (all(c("CHROM","seg_start","seg_end") %in% names(bc)))
            distinct(bc, CHROM, seg, seg_start, seg_end) else NULL

  purrr::map_dfr(seq_len(nrow(sg)), function(i) {
    r <- sg[i, ]
    tibble::tibble(
      sample_id = s, seg = r$seg, state = r[[state_col]], LLR = r$LLR %||% NA_real_,
      jp_by_name  = r$seg %in% unique(jp$seg),
      bc_by_name  = r$seg %in% unique(bc$seg),
      jp_by_coord = overlaps(r$CHROM, r$seg_start, r$seg_end, jp_u),
      bc_by_coord = overlaps(r$CHROM, r$seg_start, r$seg_end, bc_u))
  })
})

d <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
stopifnot(nrow(d) > 0)

cat("=== matching by NAME vs by COORDINATE ===\n")
print(as.data.frame(d |> summarise(
  segments            = n(),
  in_joint_post_name  = sum(jp_by_name),
  in_joint_post_coord = sum(jp_by_coord, na.rm = TRUE),
  in_bulk_name        = sum(bc_by_name),
  in_bulk_coord       = sum(bc_by_coord, na.rm = TRUE))), row.names = FALSE)

cat("\n=== segments absent by NAME but present by COORDINATE (i.e. renamed) ===\n")
ren <- d |> mutate(
  jp_renamed = !jp_by_name &  jp_by_coord,
  bc_renamed = !bc_by_name &  bc_by_coord,
  jp_truly_absent = !jp_by_name & !jp_by_coord,
  bc_truly_absent = !bc_by_name & !bc_by_coord)
print(as.data.frame(ren |> summarise(
  joint_post_renamed      = sum(jp_renamed, na.rm = TRUE),
  joint_post_truly_absent = sum(jp_truly_absent, na.rm = TRUE),
  bulk_renamed            = sum(bc_renamed, na.rm = TRUE),
  bulk_truly_absent       = sum(bc_truly_absent, na.rm = TRUE))), row.names = FALSE)

cat("\n=== the panels compared on COORDINATES rather than labels ===\n")
agree <- ren |> group_by(sample_id) |>
  summarise(mismatch_by_name  = sum(jp_by_name != bc_by_name),
            mismatch_by_coord = sum(jp_by_coord != bc_by_coord, na.rm = TRUE),
            .groups = "drop")
cat(sprintf("samples agreeing by NAME      : %d / %d\n",
            sum(agree$mismatch_by_name == 0), nrow(agree)))
cat(sprintf("samples agreeing by COORDINATE: %d / %d\n",
            sum(agree$mismatch_by_coord == 0), nrow(agree)))

cat("\n=== still disagreeing on coordinates ===\n")
bad <- agree |> filter(mismatch_by_coord > 0)
if (nrow(bad)) print(as.data.frame(bad), row.names = FALSE) else cat("  none\n")

readr::write_csv(ren, "results/segment_label_mismatch.csv")
cat("\nwrote results/segment_label_mismatch.csv\nDIAG DONE\n")
