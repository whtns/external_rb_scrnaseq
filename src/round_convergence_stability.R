# Is there a measurable convergence point in numbat's iteration?
#
# Tests whether consecutive rounds stop changing. For each adjacent pair
# (k, k+1) we compute the agreement between their non-neutral segment sets:
# a segment matches if the other round calls a same-direction segment on the
# same chromosome with >=50% reciprocal overlap. Jaccard = matched / union.
#
# If agreement jumps between r1->r2 and r2->r3 and then plateaus, "numbat has
# converged by round 2" is a defensible, data-backed statement, and round 1 is
# the pre-convergence pass rather than one option among four.
suppressPackageStartupMessages({
  library(data.table)
})

EXC    <- c("SRX11133592", "SRX11133593", "SRX11133594")
NB_DIR <- "output/numbat_sridhar"

man <- fread("results/numbat_selected_round.csv")
man <- man[grepl("^SRX", sample_id) & !sample_id %in% EXC]
samples <- man$sample_id

dir_of <- function(x) ifelse(grepl("amp|bamp", x), "gain",
                      ifelse(grepl("del|loh|bdel", x), "loss", "other"))

read_nn <- function(s, k) {
  p <- file.path(NB_DIR, s, sprintf("segs_consensus_%d.tsv", k))
  if (!file.exists(p)) return(NULL)
  d <- as.data.frame(fread(p))
  st <- if ("cnv_state_post" %in% names(d)) d$cnv_state_post else d$cnv_state
  d$state <- as.character(st)
  d$CHROM <- as.character(d$CHROM)
  d <- unique(d[, c("CHROM", "seg_start", "seg_end", "state", "LLR")])
  d <- d[!grepl("^neu", d$state), , drop = FALSE]
  if (!nrow(d)) return(d)
  d$dir <- dir_of(d$state)
  d
}

# how many rows of `a` have a >=50% reciprocal-overlap partner in `b`
n_matched <- function(a, b) {
  if (is.null(a) || !nrow(a) || is.null(b) || !nrow(b)) return(0L)
  sum(vapply(seq_len(nrow(a)), function(i) {
    cand <- b[b$CHROM == a$CHROM[i] & b$dir == a$dir[i], , drop = FALSE]
    if (!nrow(cand)) return(FALSE)
    ov <- pmax(0, pmin(cand$seg_end, a$seg_end[i]) - pmax(cand$seg_start, a$seg_start[i]))
    la <- a$seg_end[i] - a$seg_start[i]
    lb <- cand$seg_end - cand$seg_start
    any(ov / la >= 0.5 & ov / lb >= 0.5)
  }, logical(1)))
}

rows <- list()
for (s in samples) {
  ks <- as.integer(gsub("segs_consensus_|\\.tsv", "",
                        list.files(file.path(NB_DIR, s), "^segs_consensus_")))
  ks <- sort(ks[!is.na(ks)])
  if (length(ks) < 2) next
  for (j in seq_len(length(ks) - 1L)) {
    k1 <- ks[j]; k2 <- ks[j + 1L]
    a <- read_nn(s, k1); b <- read_nn(s, k2)
    if (is.null(a) || is.null(b)) next
    ma <- n_matched(a, b); mb <- n_matched(b, a)
    uni <- nrow(a) + nrow(b) - ma
    rows[[length(rows) + 1L]] <- data.table(
      sample_id = s, pair = sprintf("r%d->r%d", k1, k2), from = k1,
      n_from = nrow(a), n_to = nrow(b),
      matched_from = ma, matched_to = mb,
      jaccard = if (uni > 0) ma / uni else NA_real_,
      frac_from_kept = if (nrow(a)) ma / nrow(a) else NA_real_,
      frac_to_new    = if (nrow(b)) 1 - mb / nrow(b) else NA_real_)
  }
}
S <- rbindlist(rows)
fwrite(S, "results/round_convergence_stability.csv")

cat("=== segment-set agreement between consecutive numbat rounds ===\n")
cat("(>=50% reciprocal overlap, same chromosome, same direction)\n\n")
print(S[, .(samples = .N,
            mean_jaccard   = round(mean(jaccard, na.rm = TRUE), 3),
            median_jaccard = round(median(jaccard, na.rm = TRUE), 3),
            mean_kept      = round(mean(frac_from_kept, na.rm = TRUE), 3),
            mean_new       = round(mean(frac_to_new, na.rm = TRUE), 3)),
        keyby = pair])

cat("\npaired test, r1->r2 jaccard vs r2->r3 jaccard:\n")
a <- S[pair == "r1->r2", .(sample_id, j1 = jaccard)]
b <- S[pair == "r2->r3", .(sample_id, j2 = jaccard)]
d <- merge(a, b, by = "sample_id")
if (nrow(d) > 3) {
  print(wilcox.test(d$j1, d$j2, paired = TRUE))
  cat("  mean r1->r2 =", round(mean(d$j1, na.rm = TRUE), 3),
      "  mean r2->r3 =", round(mean(d$j2, na.rm = TRUE), 3),
      "  (n =", nrow(d), ")\n")
}
cat("\npaired test, r2->r3 jaccard vs r3->r4 jaccard:\n")
b2 <- S[pair == "r3->r4", .(sample_id, j3 = jaccard)]
d2 <- merge(b, b2, by = "sample_id")
if (nrow(d2) > 3) {
  print(wilcox.test(d2$j2, d2$j3, paired = TRUE))
  cat("  mean r2->r3 =", round(mean(d2$j2, na.rm = TRUE), 3),
      "  mean r3->r4 =", round(mean(d2$j3, na.rm = TRUE), 3),
      "  (n =", nrow(d2), ")\n")
}

# per-sample: is the round-2 segment set already the "settled" one?
cat("\nsamples whose segment set is unchanged (jaccard = 1) from that round on:\n")
print(S[, .(n_stable = sum(jaccard >= 0.999, na.rm = TRUE), n = .N), keyby = pair])
cat("\nwrote results/round_convergence_stability.csv\n")
