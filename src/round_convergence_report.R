# Evidence for a fixed-round numbat policy.
#
# Question: is there a defensible principle for choosing how many numbat
# iterations to report? This script tests three claims:
#   (1) round 1 over-calls segments relative to later rounds;
#   (2) the segments round 1 adds are disproportionately small/low-LLR,
#       i.e. the focal calls most likely to be artifactual;
#   (3) a fixed default of i=2 loses few canonical RB events, so per-sample
#       deviation can be reserved for named, justified exceptions.
#
# Canonical RB events use the SAME arm-coverage-thresholded definition as the
# reported QC table (numbat_rds_qc.csv), not the older unthresholded rule that
# results/numbat_selected_round.csv was built with.
suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

EXC      <- c("SRX11133592", "SRX11133593", "SRX11133594")
NB_DIR   <- "output/numbat_sridhar"
FOCAL_MB <- 10          # "focal" = shorter than this
MAX_K    <- 12L

man <- fread("results/numbat_selected_round.csv")
man <- man[grepl("^SRX", sample_id) & !sample_id %in% EXC]
samples <- man$sample_id
cat("samples in scope:", length(samples), "\n")

nonneu <- function(d) {
  st <- if ("cnv_state_post" %in% names(d)) d$cnv_state_post else d$cnv_state
  d$state <- as.character(st)
  d[!grepl("^neu", d$state), ]
}

## ---- per-sample x per-round summary -------------------------------------
rows <- list(); segstore <- list()
for (s in samples) {
  for (k in 1:MAX_K) {
    p <- file.path(NB_DIR, s, sprintf("segs_consensus_%d.tsv", k))
    if (!file.exists(p)) next
    d <- as.data.frame(fread(p))
    d$CHROM <- as.character(d$CHROM)
    d <- unique(d[, c("CHROM", "seg", "cnv_state", "cnv_state_post",
                      "seg_start", "seg_end", "LLR", "n_genes", "n_snps")])
    nn <- nonneu(d)
    nn$size_mb <- (nn$seg_end - nn$seg_start) / 1e6
    ev <- numbat_rb_events(d)
    segstore[[paste(s, k)]] <- data.table(sample_id = s, round = k, nn)
    rows[[length(rows) + 1L]] <- data.table(
      sample_id  = s,
      round      = k,
      n_segs     = nrow(nn),
      n_focal    = sum(nn$size_mb <  FOCAL_MB),
      n_large    = sum(nn$size_mb >= FOCAL_MB),
      median_mb  = if (nrow(nn)) median(nn$size_mb) else NA_real_,
      median_llr = if (nrow(nn)) median(nn$LLR, na.rm = TRUE) else NA_real_,
      n_rb       = length(ev),
      rb_events  = paste(sort(ev), collapse = ";"))
  }
}
byround <- rbindlist(rows)
segs    <- rbindlist(segstore)
byround[, max_round := max(round), by = sample_id]
fwrite(byround, "results/round_convergence_by_round.csv")

## ---- claim 1: does round 1 call the most segments? ----------------------
cat("\n=== segments called per round (non-neutral, deduplicated) ===\n")
print(byround[, .(samples = .N, mean_segs = round(mean(n_segs), 1),
                  median_segs = as.numeric(median(n_segs)),
                  mean_focal = round(mean(n_focal), 1),
                  median_mb = round(median(median_mb, na.rm = TRUE), 1)),
              keyby = round])

r1 <- byround[round == 1, .(sample_id, s1 = n_segs, f1 = n_focal)]
r2 <- byround[round == 2, .(sample_id, s2 = n_segs, f2 = n_focal)]
cmp <- merge(r1, r2, by = "sample_id")
cmp[, d_segs := s1 - s2][, d_focal := f1 - f2]
cat("\nround 1 vs round 2, per sample (n =", nrow(cmp), "):\n")
cat("  round 1 calls MORE segments:", sum(cmp$d_segs > 0), "samples\n")
cat("  equal                      :", sum(cmp$d_segs == 0), "samples\n")
cat("  round 1 calls FEWER        :", sum(cmp$d_segs < 0), "samples\n")
cat("  total segments  r1 =", sum(cmp$s1), " r2 =", sum(cmp$s2),
    sprintf(" (%+.1f%%)\n", 100 * (sum(cmp$s1) - sum(cmp$s2)) / sum(cmp$s2)))
cat("  focal (<", FOCAL_MB, "Mb) r1 =", sum(cmp$f1), " r2 =", sum(cmp$f2), "\n")
print(wilcox.test(cmp$s1, cmp$s2, paired = TRUE))
fwrite(cmp, "results/round1_vs_round2_segments.csv")

## ---- claim 2: are the segments round 1 drops small / low-LLR? -----------
# A round-1 segment "persists" if round 2 calls a same-direction segment on the
# same chromosome with >=50% reciprocal overlap.
dir_of <- function(x) ifelse(grepl("amp|bamp", x), "gain",
                      ifelse(grepl("del|loh|bdel", x), "loss", "other"))
pers <- list()
for (s in samples) {
  a <- segs[sample_id == s & round == 1]
  b <- segs[sample_id == s & round == 2]
  if (!nrow(a)) next
  a[, dir := dir_of(state)]
  if (nrow(b)) b[, dir := dir_of(state)]
  a[, persists := FALSE]
  for (i in seq_len(nrow(a))) {
    if (!nrow(b)) break
    cand <- b[CHROM == a$CHROM[i] & dir == a$dir[i]]
    if (!nrow(cand)) next
    ov <- pmax(0, pmin(cand$seg_end, a$seg_end[i]) - pmax(cand$seg_start, a$seg_start[i]))
    la <- a$seg_end[i] - a$seg_start[i]
    lb <- cand$seg_end - cand$seg_start
    if (any(ov / la >= 0.5 & ov / lb >= 0.5)) a$persists[i] <- TRUE
  }
  pers[[s]] <- a
}
P <- rbindlist(pers)
fwrite(P, "results/round1_segment_persistence.csv")
cat("\n=== fate of round-1 segments at round 2 ===\n")
cat("round-1 non-neutral segments:", nrow(P),
    "  persist:", sum(P$persists), sprintf("(%.0f%%)", 100 * mean(P$persists)),
    "  dropped:", sum(!P$persists), "\n\n")
print(P[, .(n = .N, median_mb = round(median(size_mb), 2),
            median_LLR = round(median(LLR, na.rm = TRUE), 1),
            median_genes = as.numeric(median(n_genes, na.rm = TRUE)),
            pct_focal = round(100 * mean(size_mb < FOCAL_MB))),
        keyby = .(persists)])
cat("\nsize (Mb), dropped vs persisting:\n")
print(wilcox.test(size_mb ~ persists, data = P))
cat("\nLLR, dropped vs persisting:\n")
print(wilcox.test(LLR ~ persists, data = P))

## ---- claim 3: what does a fixed i=2 policy cost? ------------------------
ev_at <- function(s, k) {
  r <- byround[sample_id == s & round == k]
  if (!nrow(r)) return(character(0))
  if (!nzchar(r$rb_events)) return(character(0))
  strsplit(r$rb_events, ";")[[1]]
}
pol <- list()
for (s in samples) {
  ks   <- sort(byround[sample_id == s, round])
  e2   <- ev_at(s, 2L)
  eall <- unique(unlist(lapply(ks, function(k) ev_at(s, k))))
  miss <- setdiff(eall, e2)
  # smallest round that recovers everything any round finds
  fixk <- NA_integer_
  for (k in ks) if (is.na(fixk) && all(eall %in% ev_at(s, k))) fixk <- k
  pol[[s]] <- data.table(
    sample_id      = s,
    rounds_avail   = max(ks),
    n_ev_r2        = length(e2),
    ev_r2          = paste(sort(e2), collapse = ";"),
    n_ev_union     = length(eall),
    ev_union       = paste(sort(eall), collapse = ";"),
    missed_at_r2   = paste(sort(miss), collapse = ";"),
    n_missed_at_r2 = length(miss),
    min_round_full = fixk)
}
POL <- rbindlist(pol)
POL[, needs_deviation := n_missed_at_r2 > 0]
fwrite(POL, "results/round_policy_i2_vs_union.csv")

cat("\n=== cost of a fixed i = 2 default ===\n")
cat("samples:", nrow(POL), "\n")
cat("canonical RB events at i=2      :", sum(POL$n_ev_r2), "\n")
cat("canonical RB events, best round :", sum(POL$n_ev_union), "\n")
cat("events i=2 would miss           :", sum(POL$n_missed_at_r2),
    sprintf("(%.0f%% of the union)\n",
            100 * sum(POL$n_missed_at_r2) / sum(POL$n_ev_union)))
cat("samples needing deviation       :", sum(POL$needs_deviation), "of", nrow(POL), "\n\n")
if (any(POL$needs_deviation)) {
  cat("EXCEPTIONS (a canonical RB event exists at some round but not at i=2):\n")
  print(POL[needs_deviation == TRUE,
            .(sample_id, rounds_avail, ev_r2, missed_at_r2, min_round_full)])
}
cat("\nper-event breakdown of what i=2 misses:\n")
mv <- POL[n_missed_at_r2 > 0, unlist(strsplit(missed_at_r2, ";"))]
print(sort(table(mv), decreasing = TRUE))

## a fixed i=2 vs the CURRENT per-sample selection
sel <- merge(POL, man[, .(sample_id, selected_round)], by = "sample_id")
sel[, n_ev_sel := mapply(function(s, k) length(ev_at(s, as.integer(k))),
                         sample_id, selected_round)]
cat("\nfixed i=2 vs current per-sample selection:\n")
cat("  events, fixed i=2        :", sum(sel$n_ev_r2), "\n")
cat("  events, current selection:", sum(sel$n_ev_sel), "\n")
cat("  current selection uses round 1 for", sum(sel$selected_round == 1), "of",
    nrow(sel), "samples\n")
fwrite(sel, "results/round_policy_summary.csv")
cat("\nwrote results/round_convergence_by_round.csv,",
    "round1_vs_round2_segments.csv,\n  round1_segment_persistence.csv,",
    "round_policy_i2_vs_union.csv, round_policy_summary.csv\n")
