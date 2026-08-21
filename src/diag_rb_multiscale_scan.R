# Multi-scale scan for canonical RB SCNAs within the five known loci.
# SRX samples only.
#
# WHY THIS EXISTS
#
# src/diag_rb_targeted_five_locus.R tested each canonical arm as ONE fixed
# interval and lost to the cross-round union, 83 events vs 123. The diagnosis
# was dilution, not statistics: events that test caught span a median 67% of
# their arm, events it missed span a median 36%. Averaging over a whole arm
# washes out focal sub-arm events, and much of the real RB signal is focal.
#
# So keep the prior -- only these five loci matter -- but stop assuming the event
# fills the arm. Scan nested sub-windows from 4 Mb up to the whole arm. That is
# still vastly less multiple testing than a genome-wide boundary search, and the
# within-locus search is corrected exactly (see NULL below).
#
# Reads bulk_clones_final.tsv.gz. No numbat rerun.
#
# COMMON GRID
#
# The earlier script put expression on gene indices and allelic signal on SNP
# blocks -- two different index spaces, so a single circular shift could not
# drive both and the channels had to be combined after the fact. Here both are
# binned onto the same fixed BIN_MB physical grid, so one shift moves both
# channels together and a genuinely joint statistic is available per window.
#
#   expression  mean gene-level logFC in the bin
#   allelic     |sum(pAD)/sum(DP) - 0.5| in the bin
#
# The allelic channel must be built from numbat's PHASED alt depth (pAD) at this
# bin scale, never from per-SNP BAF or per-SNP major/minor counts. Median SNP
# depth here is 2 reads, so per-SNP BAF is effectively 0/0.5/1; and per-SNP
# major/minor has a max() bias that scales with local depth, which itself rises
# with copy number, so an amplified arm reads as MORE balanced. Binning phased
# pAD avoids both. Bins are also short enough that population phase switches
# (which wash out whole-arm phased BAF) do not average the imbalance away.
#
# NULL -- MASKED to diploid territory, exact, and it corrects the window search
#
# The first version of this script shifted the window family over the WHOLE
# genome and collapsed: 16 calls, with 1q_gain going 27 -> 0 despite z_max = 6.9
# (expr_z 3.2, ai_z 6.7). The reason is that a max-over-windows null placed
# anywhere in an aneuploid tumor genome lands on that tumor's OTHER real CNVs.
# The test was effectively asking "is 1q stronger than the strongest sub-window
# anywhere else in this genome" -- a bar a real event usually fails, because RB
# tumors carry many. Contamination hurts a max statistic far more than a single
# fixed window, because every null draw gets the best of ~36 tries.
#
# So null placements are now restricted to bins numbat never called non-neutral
# in ANY consensus round for that sample. Measured across the 39 SRX samples,
# the union of non-neutral segments covers a median of just 14% of the genome
# (q75 = 22%), so ~86% remains as genuinely diploid null territory; only 1 of 39
# samples is left with under 25%. An offset counts only if the entire shifted
# locus footprint lands in clean bins.
#
# Taking the best of many windows inflates the statistic, so the null must be
# the null of that same maximum. For every circular-shift offset of the genome,
# the WHOLE window family is shifted together and the max combined z over the
# family is recomputed. The observed value is simply that max at offset zero.
# This corrects the within-locus multiplicity exactly, with no Bonferroni or
# Sidak approximation, and it automatically accounts for the heavy correlation
# between overlapping windows.
#
# Window means are standardized per window length against that length's own
# circular-shift null, so scales are comparable before the max is taken.
#
# TAIL. The empirical max-null gives p >= 1/(n_bins + 1) ~ 3e-4, too coarse to
# survive FDR over the 195-test grid. A maximum of many correlated windows is
# asymptotically Gumbel, so a Gumbel is fitted to the null max values and used
# for the tail; the empirical p is carried alongside as a floor-limited check.
# Where the two disagree badly, trust the empirical one and say so.
#
# Usage:
#   Rscript src/diag_rb_multiscale_scan.R [numbat_dir] [suffix]
#
# Output:
#   results/rb_multiscale_scan[suffix].csv

suppressPackageStartupMessages({
  library(data.table)
})

args       <- commandArgs(trailingOnly = TRUE)
NUMBAT_DIR <- if (length(args) >= 1) args[[1]] else "output/numbat_sridhar"
SUFFIX     <- if (length(args) >= 2) args[[2]] else ""

BIN_MB      <- 2        # physical bin size
MIN_GENES   <- 2        # genes needed for a bin to carry expression
MIN_SNPS    <- 15       # SNPs needed for a bin to carry allelic signal
MIN_VALID   <- 100      # clean null placements needed to report a p-value
MIN_CELLS   <- 20
SCALES_MB   <- c(4, 8, 16, 32, 64, 128)   # window lengths to scan, plus whole arm
MIN_BINS    <- 4
FDR         <- 0.05

# The five canonical RB arm events. hg38 centromere midpoints, matching
# src/diag_rb_scna_recovery.R. chr13 is acrocentric so all of it is 13q here.
LOCI <- data.table(
  event = c("1q_gain", "2p_gain", "6p_gain", "13q_loss", "16q_loss"),
  CHROM = c("1", "2", "6", "13", "16"),
  start = c(125.0e6, 0,      0,      0,   36.8e6),
  end   = c(Inf,     93.0e6, 59.0e6, Inf, Inf),
  dir   = c(1,       1,      1,      -1,  -1)      # +1 gain, -1 loss
)

CHROM_ORDER <- c(as.character(1:22), "X")

# ---------------------------------------------------------------------------
# All circular-shift window means of length k, as a vector indexed by start.
# ---------------------------------------------------------------------------
window_means <- function(v, k) {
  n  <- length(v)
  cs <- cumsum(c(0, rep(v, 2)))
  (cs[seq_len(n) + k] - cs[seq_len(n)]) / k
}

# Gumbel fit by moments: sd = beta*pi/sqrt(6), mean = mu + beta*gamma_em.
gumbel_p <- function(null_max, obs) {
  s <- stats::sd(null_max); m <- mean(null_max)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  beta <- s * sqrt(6) / pi
  mu   <- m - beta * 0.5772156649
  p    <- -expm1(-exp(-(obs - mu) / beta))   # 1 - exp(-exp(-z)), stable
  min(max(p, .Machine$double.xmin), 1)
}

# ---------------------------------------------------------------------------
# Union of every non-neutral segment across all consensus rounds for a sample.
# These regions are excluded from null placements.
# ---------------------------------------------------------------------------
read_mask <- function(dir) {
  fs <- list.files(dir, pattern = "^segs_consensus_[0-9]+\\.tsv$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  m <- rbindlist(lapply(fs, function(f) {
    z <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
    if (is.null(z) || nrow(z) == 0) return(NULL)
    st <- if ("cnv_state_post" %in% names(z)) z$cnv_state_post else z$cnv_state
    z  <- z[!is.na(st) & st != "neu"]
    if (nrow(z) == 0) return(NULL)
    data.table(CHROM = as.character(z$CHROM),
               s = as.numeric(z$seg_start), e = as.numeric(z$seg_end))
  }), fill = TRUE)
  if (is.null(m) || nrow(m) == 0) return(NULL)
  m[is.finite(s) & is.finite(e) & e > s]
}

# ---------------------------------------------------------------------------
# One sample.
# ---------------------------------------------------------------------------
test_sample <- function(sample_id, path, mask = NULL) {
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

  BIN <- BIN_MB * 1e6
  out <- list()

  for (cl in clones$sample) {
    y <- x[sample == cl]

    # --- expression per bin ------------------------------------------------
    g <- unique(y[!is.na(logFC) & !is.na(gene) & !is.na(gene_start),
                  .(gene, ci, gene_start, logFC)], by = "gene")
    g[, bin := floor(gene_start / BIN)]
    ge <- g[, .(expr = mean(logFC), ng = .N), by = .(ci, bin)][ng >= MIN_GENES]

    # --- allelic per bin ---------------------------------------------------
    s <- y[!is.na(pAD) & !is.na(DP) & DP > 0, .(ci, POS, pAD, DP)]
    s[, bin := floor(POS / BIN)]
    ab <- s[, .(ai = abs(sum(pAD) / sum(DP) - 0.5), ns = .N),
            by = .(ci, bin)][ns >= MIN_SNPS]

    # --- common grid: bins carrying BOTH channels --------------------------
    b <- merge(ge[, .(ci, bin, expr)], ab[, .(ci, bin, ai)], by = c("ci", "bin"))
    if (nrow(b) < 200L) next
    setorder(b, ci, bin)
    n <- nrow(b)

    # Clean = never called non-neutral in any consensus round. Null placements
    # are restricted to these bins so the null is diploid rather than "somewhere
    # else in this tumor's aneuploid genome".
    b[, clean := TRUE]
    if (!is.null(mask) && nrow(mask) > 0) {
      b[, `:=`(chrom = CHROM_ORDER[ci], bs = bin * BIN, be = (bin + 1) * BIN)]
      idx <- b[mask, on = .(chrom == CHROM, bs < e, be > s),
               which = TRUE, nomatch = NULL, allow.cartesian = TRUE]
      if (length(idx)) b[unique(idx), clean := FALSE]
      b[, c("chrom", "bs", "be") := NULL]
    }

    for (i in seq_len(nrow(LOCI))) {
      L   <- LOCI[i]
      lci <- match(L$CHROM, CHROM_ORDER)
      idx <- which(b$ci == lci &
                   b$bin >= floor(L$start / BIN) &
                   b$bin <= (if (is.infinite(L$end)) Inf else floor(L$end / BIN)))
      if (length(idx) < MIN_BINS) next

      a_lo  <- min(idx); a_hi <- max(idx)
      nloc  <- a_hi - a_lo + 1L

      # --- window family: scales that fit, plus the whole locus ------------
      ks <- unique(c(SCALES_MB[SCALES_MB / BIN_MB <= nloc] / BIN_MB, nloc))
      ks <- ks[ks >= MIN_BINS & ks <= floor(n / 10)]
      if (length(ks) == 0) next

      fam <- rbindlist(lapply(ks, function(k) {
        step   <- max(1L, floor(k / 2))
        starts <- seq(a_lo, a_hi - k + 1L, by = step)
        if (length(starts) == 0) return(NULL)
        data.table(start = starts, k = as.integer(k))
      }))
      if (is.null(fam) || nrow(fam) == 0) next

      # --- per-length null moments, both channels --------------------------
      Me <- Ma <- list()
      mom <- list()
      for (k in unique(fam$k)) {
        me <- window_means(b$expr, k)
        ma <- window_means(b$ai,   k)
        Me[[as.character(k)]] <- me
        Ma[[as.character(k)]] <- ma
        mom[[as.character(k)]] <- list(
          em = mean(me), es = stats::sd(me),
          am = mean(ma), as = stats::sd(ma))
      }

      # --- joint z for every window at every offset ------------------------
      # Window j at offset o sits at start (fam$start[j] + o), wrapped. Both
      # channels move together because they share the grid.
      zmax <- rep(-Inf, n)
      for (j in seq_len(nrow(fam))) {
        kk <- as.character(fam$k[j])
        mm <- mom[[kk]]
        if (!is.finite(mm$es) || mm$es <= 0 || !is.finite(mm$as) || mm$as <= 0) next
        pos <- ((fam$start[j] - 1L + seq_len(n) - 1L) %% n) + 1L
        ze  <- L$dir * (Me[[kk]][pos] - mm$em) / mm$es
        za  <- (Ma[[kk]][pos] - mm$am) / mm$as
        zmax <- pmax(zmax, (ze + za) / sqrt(2))
      }
      if (!is.finite(zmax[1])) next

      # Offset t places the locus footprint at bin (a_lo - 1 + t - 1) mod n.
      # Keep only placements landing entirely in clean (diploid) bins. The true
      # placement drops out automatically whenever the locus is itself called.
      cc    <- cumsum(c(0L, rep(as.integer(!b$clean), 2)))
      st0   <- ((a_lo - 1L) + (seq_len(n) - 1L)) %% n
      valid <- (cc[st0 + nloc + 1L] - cc[st0 + 1L]) == 0L

      obs  <- zmax[1]                  # offset 0 is the true placement
      null <- zmax[valid & is.finite(zmax)]
      if (length(null) < MIN_VALID) next
      p_emp <- (sum(null >= obs) + 1) / (length(null) + 1)
      p_gum <- gumbel_p(null, obs)

      # --- describe the winning window -------------------------------------
      best_j <- NA_integer_; best_z <- -Inf
      for (j in seq_len(nrow(fam))) {
        kk <- as.character(fam$k[j]); mm <- mom[[kk]]
        if (!is.finite(mm$es) || mm$es <= 0 || !is.finite(mm$as) || mm$as <= 0) next
        ze <- L$dir * (Me[[kk]][fam$start[j]] - mm$em) / mm$es
        za <- (Ma[[kk]][fam$start[j]] - mm$am) / mm$as
        zc <- (ze + za) / sqrt(2)
        if (zc > best_z) { best_z <- zc; best_j <- j }
      }
      kk <- as.character(fam$k[best_j]); mm <- mom[[kk]]
      w_ze <- L$dir * (Me[[kk]][fam$start[best_j]] - mm$em) / mm$es
      w_za <- (Ma[[kk]][fam$start[best_j]] - mm$am) / mm$as
      wb   <- b[fam$start[best_j]:(fam$start[best_j] + fam$k[best_j] - 1L)]

      # Bins are sparse (only those carrying both channels survive), so a
      # window of k bins does NOT span k * BIN_MB. Report the physical span.
      loc_bins <- b$bin[idx]
      win_mb_phys <- (max(wb$bin) + 1L - min(wb$bin)) * BIN / 1e6
      loc_mb_phys <- (max(loc_bins) + 1L - min(loc_bins)) * BIN / 1e6

      out[[length(out) + 1L]] <- data.table(
        sample_id  = sample_id,
        event      = L$event,
        clone      = as.character(cl),
        n_cells    = clones[sample == cl, n_cells][1],
        n_bins_loc = nloc,
        n_windows  = nrow(fam),
        n_null     = length(null),
        win_mb     = win_mb_phys,
        win_start  = min(wb$bin) * BIN / 1e6,
        win_end    = (max(wb$bin) + 1L) * BIN / 1e6,
        arm_frac   = win_mb_phys / loc_mb_phys,
        expr_z     = w_ze,
        ai_z       = w_za,
        z_max      = obs,
        p_emp      = p_emp,
        p_gumbel   = p_gum
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
  per[[j]] <- test_sample(sid, file.path(dirs[j], "bulk_clones_final.tsv.gz"),
                          mask = read_mask(dirs[j]))
  invisible(gc(verbose = FALSE))
}
res <- rbindlist(per[!vapply(per, is.null, logical(1))])
if (nrow(res) == 0) stop("no testable sample x locus x clone combinations")

# --- best clone per sample x locus, Sidak for the clone search --------------
res[, n_clones_tested := .N, by = .(sample_id, event)]
setorder(res, sample_id, event, p_gumbel, na.last = TRUE)
best <- res[, .SD[1], by = .(sample_id, event)]
best[, p_sidak := 1 - (1 - p_gumbel)^pmax(n_clones_tested, 1)]

# --- BH over the sample x locus grid ---------------------------------------
best[, q := p.adjust(p_sidak, method = "BH")]
best[, called := !is.na(q) & q < FDR]

# Gains must raise expression; losses get no veto, because copy-neutral LOH
# gives allelic imbalance with flat expression. See the sibling script.
best[, is_gain := grepl("_gain$", event)]
best[, expr_consistent := is.na(expr_z) | expr_z > 0]
best[, called := called & (!is_gain | expr_consistent)]

setorder(best, sample_id, event)
dir.create("results", showWarnings = FALSE)
outfile <- file.path("results", paste0("rb_multiscale_scan", SUFFIX, ".csv"))
fwrite(best, outfile)

# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------
cat("\n=== multi-scale scan within the five canonical loci, SRX samples ===\n")
cat("samples:", uniqueN(best$sample_id), "  sample x locus tests:", nrow(best),
    "  bin:", BIN_MB, "Mb  FDR:", FDR, "\n")
cat("windows scanned per locus (median):", stats::median(best$n_windows), "\n\n")

cat("calls by event:\n")
print(best[, .(tested = .N, called = sum(called)), by = event][order(event)])

cat("\nwinning window size among calls (fraction of the arm):\n")
print(best[called == TRUE, .(n = .N,
       median_win_mb = round(stats::median(win_mb), 1),
       median_arm_frac = round(stats::median(arm_frac), 2)), by = event][order(event)])

cat("\ntotal canonical events called (multi-scale):", sum(best$called), "\n")

cmp <- "results/rb_scna_canonical_by_sample.csv"
if (file.exists(cmp)) {
  seg <- fread(cmp)
  seg <- seg[grepl("^SRX", sample_id) & mixed_provenance == FALSE &
             sample_id %in% best$sample_id]
  seg[, seg_union := called_final == TRUE | recovered_by_union == TRUE]
  cat("\n--- same samples, for comparison ---\n")
  cat("  numbat final round :", seg[called_final == TRUE, .N], "\n")
  cat("  cross-round union  :", seg[seg_union == TRUE, .N], "\n")
  cat("  whole-arm targeted :",
      if (file.exists("results/rb_targeted_five_locus.csv"))
        fread("results/rb_targeted_five_locus.csv")[called == TRUE, .N] else NA, "\n")
  cat("  multi-scale (this) :", sum(best$called), "\n")

  m <- merge(best[, .(sample_id, event, called)],
             seg[, .(sample_id, event, seg_union)],
             by = c("sample_id", "event"), all = TRUE)
  m[is.na(called), called := FALSE]
  m[is.na(seg_union), seg_union := FALSE]
  cat("\n  agreement grid:\n")
  print(table(multiscale = m$called, union = m$seg_union))
  cat("\n  multiscale-only (union misses these):", m[called & !seg_union, .N], "\n")
  cat("  union-only (multiscale misses these):", m[!called & seg_union, .N], "\n")
  cat("  combined multiscale OR union        :", m[called | seg_union, .N], "\n")
}

cat("\nwrote", outfile, "\n")
