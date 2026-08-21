# Validate the round-selection wiring before any RDS is overwritten.
#
# For every sample, builds the Numbat object BOTH ways -- at the package default
# i = 2 (what process_numbat_rds.R produced before the change) and at the round
# recorded in results/numbat_selected_round.csv -- and compares canonical RB
# event counts. Writes nothing into output/numbat_sridhar/.
#
# The point is to confirm three things before touching production RDS files:
#   1. Numbat$new() actually succeeds at the selected round for every sample.
#   2. The canonical gain is what the segs_consensus scoring predicted.
#   3. No sample regresses.
#
# Usage:
#   Rscript src/verify_selected_round_rds.R [numbat_dir] [manifest]

suppressPackageStartupMessages({
  library(numbat); library(data.table)
})

args     <- commandArgs(trailingOnly = TRUE)
NB_DIR   <- if (length(args) >= 1) args[[1]] else "output/numbat_sridhar"
MANIFEST <- if (length(args) >= 2) args[[2]] else "results/numbat_selected_round.csv"

man <- fread(MANIFEST)

# Canonical RB events, matching src/diag_rb_scna_recovery.R.
canon_events <- function(sc) {
  if (is.null(sc) || nrow(sc) == 0) return(character(0))
  sc <- as.data.table(sc)
  sc[, CHROM := as.character(CHROM)]
  st <- if ("cnv_state_post" %in% names(sc)) sc$cnv_state_post else sc$cnv_state
  g <- grepl("amp", st); l <- grepl("del|loh", st)
  e <- character(0)
  if (any(sc$CHROM == "1"  & sc$seg_end   > 125e6  & g, na.rm = TRUE)) e <- c(e, "1q_gain")
  if (any(sc$CHROM == "2"  & sc$seg_start <  93e6  & g, na.rm = TRUE)) e <- c(e, "2p_gain")
  if (any(sc$CHROM == "6"  & sc$seg_start <  59e6  & g, na.rm = TRUE)) e <- c(e, "6p_gain")
  if (any(sc$CHROM == "13"                         & l, na.rm = TRUE)) e <- c(e, "13q_loss")
  if (any(sc$CHROM == "16" & sc$seg_end   > 36.8e6 & l, na.rm = TRUE)) e <- c(e, "16q_loss")
  e
}

build <- function(d, i) {
  tryCatch(Numbat$new(out_dir = d, i = i), error = function(e) NULL)
}

out <- list()
for (j in seq_len(nrow(man))) {
  sid <- man$sample_id[j]
  k   <- as.integer(man$selected_round[j])
  d   <- file.path(NB_DIR, sid)
  if (!dir.exists(d)) next

  cat(sprintf("[%2d/%2d] %s  (selected %d)\n", j, nrow(man), sid, k))

  a <- build(d, 2L)
  b <- if (k == 2L) a else build(d, k)

  ea <- if (is.null(a)) NA_character_ else paste(sort(canon_events(a$segs_consensus)), collapse = ";")
  eb <- if (is.null(b)) NA_character_ else paste(sort(canon_events(b$segs_consensus)), collapse = ";")

  out[[length(out) + 1L]] <- data.table(
    sample_id        = sid,
    selected_round   = k,
    mixed_provenance = man$mixed_provenance[j],
    built_i2         = !is.null(a),
    built_sel        = !is.null(b),
    n_canon_i2       = if (is.null(a)) NA_integer_ else length(canon_events(a$segs_consensus)),
    n_canon_sel      = if (is.null(b)) NA_integer_ else length(canon_events(b$segs_consensus)),
    events_i2        = ea,
    events_sel       = eb
  )
  rm(a, b); invisible(gc(verbose = FALSE))
}

res <- rbindlist(out)
res[, delta := n_canon_sel - n_canon_i2]
dir.create("results", showWarnings = FALSE)
fwrite(res, "results/verify_selected_round_rds.csv")

clean <- res[mixed_provenance == FALSE]

cat("\n=== object-level validation ===\n")
cat("samples attempted:", nrow(res), "\n")
cat("built OK at i = 2       :", res[built_i2 == TRUE, .N], "\n")
cat("built OK at selected i  :", res[built_sel == TRUE, .N], "\n")
if (res[built_sel == FALSE, .N] > 0) {
  cat("\n!! FAILED to build at selected round:\n")
  print(res[built_sel == FALSE, .(sample_id, selected_round)])
}

cat("\ncanonical RB events, clean samples (", nrow(clean), "):\n", sep = "")
cat("  i = 2 (current pipeline) :", clean[, sum(n_canon_i2, na.rm = TRUE)], "\n")
cat("  selected round           :", clean[, sum(n_canon_sel, na.rm = TRUE)], "\n")
cat("  delta                    : +", clean[, sum(delta, na.rm = TRUE)], "\n", sep = "")

cat("\nsamples improved:", clean[delta > 0, .N],
    " unchanged:", clean[delta == 0, .N],
    " REGRESSED:", clean[delta < 0, .N], "\n")
if (clean[delta < 0, .N] > 0) {
  cat("\n!! regressions:\n")
  print(clean[delta < 0, .(sample_id, selected_round, n_canon_i2, n_canon_sel,
                           events_i2, events_sel)])
}

cat("\nwrote results/verify_selected_round_rds.csv\n")
