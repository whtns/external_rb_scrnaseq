# Characterise the canonical RB events that a fixed i = 2 default would miss.
#
# The distinction that matters for the manuscript: at round 2, is the event
#   (a) ABSENT   - no same-direction call anywhere on the arm; or
#   (b) SUB-THRESHOLD - the arm IS called, but the supporting segments cover
#       less than RB_MIN_ARM_FRAC of it (i.e. fragmented, not missing).
# (b) is a resolution difference between rounds, not a disagreement about
# whether the arm is altered, and is far easier to defend as a named deviation.
suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

NB_DIR <- "output/numbat_sridhar"
POL    <- fread("results/round_policy_i2_vs_union.csv")

read_round <- function(s, k) {
  p <- file.path(NB_DIR, s, sprintf("segs_consensus_%d.tsv", k))
  if (!file.exists(p)) return(NULL)
  d <- as.data.frame(fread(p))
  d$CHROM <- as.character(d$CHROM)
  unique(d[, c("CHROM", "seg", "cnv_state", "cnv_state_post",
               "seg_start", "seg_end", "LLR", "n_genes", "n_snps")])
}

# max LLR among the segments supporting one canonical event, at one round
supp <- function(d, event) {
  w <- .RB_WINDOWS[.RB_WINDOWS$event == event, ]
  if (is.null(d) || !nrow(d)) return(list(llr = NA_real_, n = 0L, mb = 0))
  st <- if ("cnv_state_post" %in% names(d)) d$cnv_state_post else d$cnv_state
  keep <- d$CHROM == w$CHROM &
    grepl(if (w$want == "gain") "amp" else "del|loh", as.character(st)) &
    d$seg_end > w$xmin * 1e6 & d$seg_start < w$xmax * 1e6
  if (!any(keep, na.rm = TRUE)) return(list(llr = NA_real_, n = 0L, mb = 0))
  s <- d[which(keep), , drop = FALSE]
  list(llr = max(s$LLR, na.rm = TRUE), n = nrow(s),
       mb = sum(pmin(s$seg_end, w$xmax * 1e6) - pmax(s$seg_start, w$xmin * 1e6)) / 1e6)
}

rows <- list()
for (i in seq_len(nrow(POL))) {
  s <- POL$sample_id[i]
  if (!nzchar(POL$missed_at_r2[i])) next
  evs <- strsplit(POL$missed_at_r2[i], ";")[[1]]
  d2  <- read_round(s, 2L)
  af2 <- if (is.null(d2)) NULL else numbat_rb_arm_frac(d2)
  # gsub, not sub: the filename has BOTH the prefix and the .tsv suffix to strip,
  # and sub() would drop only the first, leaving "1.tsv" -> NA.
  ks  <- as.integer(gsub("segs_consensus_|\\.tsv", "",
                         list.files(file.path(NB_DIR, s), "^segs_consensus_")))
  ks  <- sort(ks[!is.na(ks)])
  for (e in evs) {
    # first round that calls this event above threshold
    kk <- NA_integer_; afk <- NA_real_; sk <- NULL
    for (k in ks) {
      dk <- read_round(s, k)
      if (is.null(dk)) next
      if (e %in% numbat_rb_events(dk)) {
        kk <- k; afk <- numbat_rb_arm_frac(dk)[[e]]; sk <- supp(dk, e); break
      }
    }
    s2 <- supp(d2, e)
    rows[[length(rows) + 1L]] <- data.table(
      sample_id     = s,
      event         = e,
      arm_frac_r2   = if (is.null(af2)) NA_real_ else round(af2[[e]], 3),
      n_supp_r2     = s2$n,
      mb_supp_r2    = round(s2$mb, 2),
      max_LLR_r2    = round(s2$llr, 1),
      recover_round = kk,
      arm_frac_rk   = round(afk, 3),
      n_supp_rk     = if (is.null(sk)) NA_integer_ else sk$n,
      mb_supp_rk    = if (is.null(sk)) NA_real_ else round(sk$mb, 2),
      max_LLR_rk    = if (is.null(sk)) NA_real_ else round(sk$llr, 1))
  }
}
M <- rbindlist(rows)
M[, status := ifelse(is.na(arm_frac_r2) | arm_frac_r2 == 0, "absent at r2",
                     "sub-threshold at r2")]
setorder(M, status, -arm_frac_r2)
fwrite(M, "results/round_policy_missed_events.csv")

cat("=== the", nrow(M), "canonical RB events a fixed i=2 would miss ===\n\n")
print(M[, .(n = .N,
            median_arm_frac_r2 = round(median(arm_frac_r2, na.rm = TRUE), 3),
            median_LLR_r2 = round(median(max_LLR_r2, na.rm = TRUE), 1),
            median_arm_frac_rk = round(median(arm_frac_rk), 3),
            median_LLR_rk = round(median(max_LLR_rk), 1)),
        keyby = status])
cat("\nby event:\n"); print(dcast(M, event ~ status, fun.aggregate = length))
cat("\nrecovering round:\n"); print(table(M$recover_round))
cat("\nfull table:\n"); print(M)
cat("\nwrote results/round_policy_missed_events.csv\n")
