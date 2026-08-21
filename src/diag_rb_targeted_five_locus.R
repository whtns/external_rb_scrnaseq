# Targeted five-locus test for canonical RB SCNAs. SRX samples only.
#
# WHY THIS EXISTS
#
# numbat searches the genome for segment boundaries and then asks whether each
# discovered segment is non-neutral. That throws away a strong prior: in
# retinoblastoma the events of interest are five known arm-level intervals. A
# genome-wide search pays the multiple-testing cost of every candidate boundary;
# testing five fixed intervals pays for five hypotheses, so the same evidence
# clears a much lower bar.
#
# Measured motivation (src/diag_rb_scna_recovery.R, 33 clean 4-round samples):
# canonical events called per consensus round run 106 / 98 / 92 / 86 while the
# canonical FRACTION of all calls stays pinned at ~26%. No numbat parameter is
# selective for these events, because numbat has no notion of which events are
# canonical. So stop asking it to discover them and test for them directly.
#
# Reads bulk_clones_final.tsv.gz, which numbat already wrote for every sample.
# No rerun. "final" is the LAST consensus round, i.e. the round that calls the
# FEWEST canonical events -- testing the worst round is the conservative choice.
#
# THE TEST -- two channels per clone per locus
#
# expression  mean gene-level logFC over the interval, signed by event direction
#             (gains expect positive, losses negative).
#
# allelic     phase-corrected allelic imbalance. Per-SNP depth in this cohort is
#             tiny (median DP = 2), so per-SNP BAF is effectively 0/0.5/1 and
#             carries no usable signal on its own; and taking per-SNP
#             major/minor counts is worse than useless, because the max() bias
#             depends on local depth, which itself rises with copy number. So we
#             aggregate numbat's PHASED alt depth (pAD) over blocks of
#             BLOCK_SNPS consecutive SNPs and take |sum(pAD)/sum(DP) - 0.5| per
#             block. Blocking matters: population phasing switches haplotype
#             over long spans, so summing pAD across a whole arm averages the
#             imbalance back to 0.5 (measured: |dev| = 0.02 whole-arm vs 0.099
#             blocked on a known 6p gain). Imbalance is direction-agnostic --
#             amp, del and loh all raise it -- so this channel corroborates and
#             the expression channel assigns direction.
#
# NULL for both channels: CIRCULAR SHIFT along the genome-ordered vector. Genes
# on one arm are strongly correlated -- that correlation is the signal, not a
# nuisance -- so a t-test or Wilcoxon against "all other genes" is badly
# anticonservative. A circular shift preserves the autocorrelation exactly and
# asks how unusual a contiguous window of this length is anywhere in the genome.
# It is also self-calibrating per clone, which matters because the noise level
# varies several-fold across clones of different size.
#
# CAVEAT ON p-VALUES. The exact permutation p is floored at 1/(n_offsets + 1),
# which for the blocked allelic channel is ~1e-3 -- too coarse to survive FDR
# over the sample x locus grid. So the reported q-values are built from the
# standardized z against the empirical null (normal tail), with the exact
# permutation p carried alongside as a sanity check. The blocked allelic null is
# right-skewed (bounded at 0), so its normal-tail p is mildly anticonservative
# in the extreme tail; treat a lone allelic call near the threshold with care.
# The two channels are combined by Stouffer, which requires both to point the
# right way to accumulate evidence.
#
# Usage:
#   Rscript src/diag_rb_targeted_five_locus.R [numbat_dir] [suffix]
#
# Output:
#   results/rb_targeted_five_locus[suffix].csv

suppressPackageStartupMessages({
  library(data.table)
})

args       <- commandArgs(trailingOnly = TRUE)
NUMBAT_DIR <- if (length(args) >= 1) args[[1]] else "output/numbat_sridhar"
SUFFIX     <- if (length(args) >= 2) args[[2]] else ""

MIN_CELLS   <- 20    # clones smaller than this are too noisy to test
MIN_GENES   <- 40    # measured genes needed in the interval
MIN_BLOCKS  <- 8     # allelic blocks needed in the interval
BLOCK_SNPS  <- 100   # SNPs per phasing-robust block
FDR         <- 0.05

# The five canonical RB arm events. hg38 centromere midpoints, matching the
# boundaries documented in src/diag_rb_scna_recovery.R. chr13 is acrocentric so
# the whole chromosome is 13q for this purpose.
LOCI <- data.table(
  event = c("1q_gain", "2p_gain", "6p_gain", "13q_loss", "16q_loss"),
  CHROM = c("1", "2", "6", "13", "16"),
  start = c(125.0e6, 0,      0,      0,   36.8e6),
  end   = c(Inf,     93.0e6, 59.0e6, Inf, Inf),
  dir   = c(1,       1,      1,      -1,  -1)      # +1 gain, -1 loss
)

CHROM_ORDER <- c(as.character(1:22), "X")

# ---------------------------------------------------------------------------
# Circular-shift null over a genome-ordered vector.
#
# v     genome-ordered values, no NAs
# k     window length (the interval occupies k contiguous positions)
# obs   observed window mean
# side  +1 upper tail, -1 lower tail
#
# All n offsets cost O(n) via a doubled cumulative sum.
# ---------------------------------------------------------------------------
circshift_null <- function(v, k, obs, side = 1) {
  n <- length(v)
  if (k < 2L || n < 10L * k) return(NULL)

  cs   <- cumsum(c(0, rep(v, 2)))
  s    <- seq_len(n)
  null <- (cs[s + k] - cs[s]) / k
  null <- null[is.finite(null)]
  if (length(null) < 50L) return(NULL)

  mu <- mean(null)
  sd <- stats::sd(null)
  z  <- if (is.finite(sd) && sd > 0) (obs - mu) / sd else NA_real_

  # Exact permutation p, floored at 1/(n+1); see the caveat in the header.
  p_emp <- if (side > 0) (sum(null >= obs) + 1) / (length(null) + 1)
           else          (sum(null <= obs) + 1) / (length(null) + 1)

  list(mu = mu, sd = sd, z = z, p_emp = p_emp, n_null = length(null))
}

# ---------------------------------------------------------------------------
# One sample: every clone x locus.
# ---------------------------------------------------------------------------
test_sample <- function(sample_id, path) {
  x <- tryCatch(
    fread(path, select = c("CHROM", "POS", "gene", "gene_start", "logFC",
                           "pAD", "DP", "sample", "n_cells"),
          showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)

  x[, CHROM := as.character(CHROM)]
  x <- x[CHROM %in% CHROM_ORDER]
  if (nrow(x) == 0) return(NULL)
  x[, ci := match(CHROM, CHROM_ORDER)]

  clones <- unique(x[, .(sample, n_cells)])[n_cells >= MIN_CELLS]
  if (nrow(clones) == 0) return(NULL)

  out <- list()

  for (cl in clones$sample) {
    y <- x[sample == cl]

    # --- expression: one row per gene, genome-ordered -----------------------
    g <- unique(y[!is.na(logFC) & !is.na(gene) & !is.na(gene_start),
                  .(gene, ci, gene_start, logFC)], by = "gene")
    setorder(g, ci, gene_start)

    # --- allelic: blocks of consecutive SNPs, genome-ordered ----------------
    s <- y[!is.na(pAD) & !is.na(DP) & DP > 0, .(ci, POS, pAD, DP)]
    setorder(s, ci, POS)
    s[, blk := (seq_len(.N) - 1L) %/% BLOCK_SNPS, by = ci]
    b <- s[, .(pos_mid = as.numeric(stats::median(POS)),
               dev     = abs(sum(pAD) / sum(DP) - 0.5)),
           by = .(ci, blk)]
    setorder(b, ci, blk)

    for (i in seq_len(nrow(LOCI))) {
      L   <- LOCI[i]
      lci <- match(L$CHROM, CHROM_ORDER)

      gi <- which(g$ci == lci & g$gene_start >= L$start & g$gene_start <= L$end)
      bi <- which(b$ci == lci & b$pos_mid    >= L$start & b$pos_mid    <= L$end)

      e <- NULL
      if (length(gi) >= MIN_GENES) {
        obs <- mean(g$logFC[gi])
        e   <- circshift_null(g$logFC, length(gi), obs, side = L$dir)
        if (!is.null(e)) e$obs <- obs
      }

      a <- NULL
      if (length(bi) >= MIN_BLOCKS) {
        obs <- mean(b$dev[bi])
        a   <- circshift_null(b$dev, length(bi), obs, side = 1)
        if (!is.null(a)) a$obs <- obs
      }

      # Signed by event direction: positive means "consistent with this event".
      expr_z_dir <- if (is.null(e)) NA_real_ else L$dir * e$z

      out[[length(out) + 1L]] <- data.table(
        sample_id  = sample_id,
        event      = L$event,
        clone      = as.character(cl),
        n_cells    = clones[sample == cl, n_cells][1],
        n_genes    = length(gi),
        n_blocks   = length(bi),
        expr_mean  = if (is.null(e)) NA_real_ else e$obs,
        expr_null  = if (is.null(e)) NA_real_ else e$mu,
        expr_z     = expr_z_dir,
        expr_p_emp = if (is.null(e)) NA_real_ else e$p_emp,
        ai_mean    = if (is.null(a)) NA_real_ else a$obs,
        ai_null    = if (is.null(a)) NA_real_ else a$mu,
        ai_z       = if (is.null(a)) NA_real_ else a$z,
        ai_p_emp   = if (is.null(a)) NA_real_ else a$p_emp
      )
    }
  }
  if (length(out) == 0) return(NULL)
  rbindlist(out)
}

# ---------------------------------------------------------------------------
# Run over SRX samples.
# ---------------------------------------------------------------------------
dirs <- list.dirs(NUMBAT_DIR, recursive = FALSE)
dirs <- dirs[grepl("/SRX[0-9]+$", dirs)]
dirs <- dirs[file.exists(file.path(dirs, "bulk_clones_final.tsv.gz"))]

cat("SRX samples with bulk_clones_final.tsv.gz:", length(dirs), "\n")
if (length(dirs) == 0) stop("nothing to do")

per <- vector("list", length(dirs))
for (j in seq_along(dirs)) {
  sid <- basename(dirs[j])
  cat(sprintf("  [%2d/%2d] %s\n", j, length(dirs), sid))
  per[[j]] <- test_sample(sid, file.path(dirs[j], "bulk_clones_final.tsv.gz"))
  # These pseudobulk tables are ~440k rows each; drop them between samples so
  # peak RSS stays flat rather than growing across the cohort.
  invisible(gc(verbose = FALSE))
}
res <- rbindlist(per[!vapply(per, is.null, logical(1))])
if (nrow(res) == 0) stop("no testable sample x locus x clone combinations")

# --- Stouffer over the two channels ----------------------------------------
# Both z are oriented so that positive supports the event. If one channel is
# unavailable the other stands alone.
res[, z_comb := fifelse(
  !is.na(expr_z) & !is.na(ai_z), (expr_z + ai_z) / sqrt(2),
  fifelse(!is.na(expr_z), expr_z, ai_z))]
res[, p_comb := stats::pnorm(z_comb, lower.tail = FALSE)]

# --- best clone per sample x locus, Sidak-corrected for the clone search ----
res[, n_clones_tested := sum(!is.na(z_comb)), by = .(sample_id, event)]
setorder(res, sample_id, event, p_comb, na.last = TRUE)
best <- res[, .SD[1], by = .(sample_id, event)]
best[, p_comb_sidak := 1 - (1 - p_comb)^pmax(n_clones_tested, 1)]

# --- BH across the sample x locus grid --------------------------------------
best[, q_comb := p.adjust(p_comb_sidak, method = "BH")]
best[, called := !is.na(q_comb) & q_comb < FDR]

# Direction consistency, GAINS ONLY. An amplification must raise expression, so
# a "gain" whose expression channel points the other way is not a gain.
#
# Losses get no such veto, because copy-neutral LOH produces strong allelic
# imbalance and NO expression change; vetoing on expression would discard those
# events by construction. Measured: with the veto applied to losses, 13q_loss
# called 0/39 despite SRX11133585 showing ai_z = 7.3 at q = 1e-4.
#
# Note what this does NOT imply. Biallelic RB1 inactivation does not require 13q
# copy loss -- two point mutations, or mutation plus promoter hypermethylation,
# inactivate RB1 with no copy change and no LOH, and MYCN-amplified RB is RB1
# proficient. So a low 13q_loss rate is not by itself evidence of a broken test;
# some of these tumors have no 13q copy event to find. The veto fix is justified
# only where imbalance IS present and expression is flat.
#
# The expr_consistent column below records the direction either way, so a
# stricter downstream filter remains possible.
best[, is_gain := grepl("_gain$", event)]
best[, expr_consistent := is.na(expr_z) | expr_z > 0]
best[, called := called & (!is_gain | expr_consistent)]

# Which channel is carrying the call?
best[, support := fifelse(!called, "",
  fifelse(!is.na(expr_z) & expr_z > 1.96 & !is.na(ai_z) & ai_z > 1.96, "both",
  fifelse(!is.na(expr_z) & expr_z > 1.96, "expr",
  fifelse(!is.na(ai_z)   & ai_z   > 1.96, "allelic", "joint"))))]

setorder(best, sample_id, event)
dir.create("results", showWarnings = FALSE)
outfile <- file.path("results", paste0("rb_targeted_five_locus", SUFFIX, ".csv"))
fwrite(best, outfile)

# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------
cat("\n=== targeted five-locus test, SRX samples ===\n")
cat("samples:", uniqueN(best$sample_id),
    "  sample x locus tests:", nrow(best),
    "  FDR:", FDR, "\n\n")

cat("calls by event:\n")
print(best[, .(tested = .N, called = sum(called)), by = event][order(event)])

cat("\nchannel carrying each call:\n")
print(best[called == TRUE, .N, by = support][order(-N)])

cat("\ntotal canonical events called (targeted):", sum(best$called), "\n")

cmp <- "results/rb_scna_canonical_by_sample.csv"
if (file.exists(cmp)) {
  seg <- fread(cmp)
  seg <- seg[grepl("^SRX", sample_id) & mixed_provenance == FALSE &
             sample_id %in% best$sample_id]
  seg[, seg_union := called_final == TRUE | recovered_by_union == TRUE]

  cat("\n--- same samples, segmentation-based (src/diag_rb_scna_recovery.R) ---\n")
  cat("  numbat final round :", seg[called_final == TRUE, .N], "\n")
  cat("  cross-round union  :", seg[seg_union == TRUE, .N], "\n")
  cat("  targeted (this)    :", sum(best$called), "\n")

  m <- merge(best[, .(sample_id, event, called)],
             seg[, .(sample_id, event, seg_union)],
             by = c("sample_id", "event"), all = TRUE)
  m[is.na(called), called := FALSE]
  m[is.na(seg_union), seg_union := FALSE]
  cat("\n  agreement grid:\n")
  print(table(targeted = m$called, union = m$seg_union))
  cat("\n  targeted-only (union misses these):", m[called & !seg_union, .N], "\n")
  cat("  union-only (targeted misses these):", m[!called & seg_union, .N], "\n")
}

cat("\nwrote", outfile, "\n")
