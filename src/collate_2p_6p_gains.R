# Which samples legitimately carry 2p gain or 6p gain, post-revamp.
#
# "Legitimate" here is the majority-of-rounds rule now in force (see
# src/select_round_by_majority.R and docs/numbat_round_policy.md), not the
# old event-maximising union:
#
#   an event is REPORTED when a majority of a sample's available consensus
#   rounds call it, where a round calls it only if the supporting amp/bamp
#   segments cover >= RB_MIN_ARM_FRAC (0.15) of the target arm.
#
# The arm floor is what separates an arm-level gain from a focal speck; the
# majority floor is what stops a single unconverged round from deciding.
# Both must hold. This script adds the evidence behind each call so a reader
# can see WHY it passed: arm coverage, union extent, max LLR, phi_mle, and the
# fraction of cells posterior-positive for the gain.
#
# Reads only. Nothing writes into output/numbat_sridhar/.
#
# Outputs:
#   results/rb_2p_6p_gain_by_sample.csv   every sample x {2p,6p}, called or not
#   results/rb_2p_6p_gain_legitimate.csv  the called rows only
#
# Usage: Rscript src/collate_2p_6p_gains.R

suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

NB_DIR <- "output/numbat_sridhar"

# Same windows numbat_rb_arm_frac() uses, so coverage and segment evidence
# cannot disagree about which segments support a call.
WIN <- data.table(
  event = c("2p_gain", "6p_gain"),
  CHROM = c("2", "6"),
  arm   = c("2p", "6p"),
  lo    = c(0, 0) * 1e6,
  hi    = c(93.9, 59.8) * 1e6
)

MAJ <- fread("results/round_matching_majority.csv", colClasses = list(character = "majority_set"))
BY  <- fread("results/round_convergence_by_round.csv")

union_len <- function(lo, hi) {
  if (!length(lo)) return(0)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  tot <- 0; cs <- lo[1]; ce <- hi[1]
  for (i in seq_along(lo)[-1]) {
    if (lo[i] > ce) { tot <- tot + (ce - cs); cs <- lo[i]; ce <- hi[i] }
    else ce <- max(ce, hi[i])
  }
  tot + (ce - cs)
}

state_of <- function(d) {
  if ("cnv_state_post" %in% names(d)) {
    fifelse(is.na(d$cnv_state_post) | d$cnv_state_post == "",
            as.character(d$cnv_state), as.character(d$cnv_state_post))
  } else as.character(d$cnv_state)
}

rows <- list()
for (i in seq_len(nrow(MAJ))) {
  s  <- MAJ$sample_id[i]
  k  <- MAJ$best_round[i]
  maj <- if (nzchar(MAJ$majority_set[i])) strsplit(MAJ$majority_set[i], ";")[[1]] else character(0)

  segf <- file.path(NB_DIR, s, sprintf("segs_consensus_%d.tsv", k))
  if (!file.exists(segf)) { warning("missing ", segf); next }
  segs <- fread(segf, colClasses = list(character = "CHROM"), showProgress = FALSE)
  frac <- numbat_rb_arm_frac(segs)
  segs[, state := state_of(segs)]

  # per-cell posterior, for the cell fraction carrying the gain
  jp <- NULL
  jpf <- file.path(NB_DIR, s, sprintf("joint_post_%d.tsv", k))
  if (file.exists(jpf)) {
    jp <- tryCatch(fread(jpf, colClasses = list(character = "CHROM"), showProgress = FALSE,
                         select = c("cell", "CHROM", "seg", "seg_start", "seg_end",
                                    "p_amp", "p_bamp")),
                   error = function(e) NULL)
  }

  for (w in seq_len(nrow(WIN))) {
    ev <- WIN$event[w]
    sup <- segs[CHROM == WIN$CHROM[w] & grepl("amp", state) &
                  seg_end > WIN$lo[w] & seg_start < WIN$hi[w]]

    # Join joint_post to the supporting segments by COORDINATE OVERLAP, not by
    # seg name: the consensus re-labels segments, so segs_consensus "2b" is
    # joint_post "2a" in SRX14116947 and 6b/6d/6f is 6a/6c/6e in SRX11133588.
    # Matching on the name silently returns nothing.
    #
    # p_amp alone is also not enough. A balanced amplification is called bamp
    # and carries its posterior in p_bamp; scoring bamp segments on p_amp
    # reports 0% of cells for a clonal gain (SRX10264519/20, SRX14116946).
    cellfrac <- NA_real_
    if (!is.null(jp) && nrow(sup)) {
      j <- jp[CHROM == WIN$CHROM[w]]
      if (nrow(j)) {
        hit <- rep(FALSE, nrow(j))
        for (r in seq_len(nrow(sup))) {
          hit <- hit | (j$seg_end > sup$seg_start[r] & j$seg_start < sup$seg_end[r])
        }
        j <- j[hit]
        if (nrow(j)) {
          j[, p_gain := pmax(fcoalesce(p_amp, 0), fcoalesce(p_bamp, 0))]
          cellfrac <- j[, .(pos = any(p_gain > 0.5, na.rm = TRUE)), by = cell][, mean(pos)]
        }
      }
    }

    nr  <- BY[sample_id == s, .N]
    ncall <- BY[sample_id == s, sum(vapply(rb_events, function(z)
      ev %in% strsplit(z, ";")[[1]], logical(1)))]

    rows[[length(rows) + 1L]] <- data.table(
      sample_id      = s,
      event          = ev,
      called         = ev %in% maj,
      selected_round = k,
      arm_frac       = round(unname(frac[ev]), 3),
      n_seg          = nrow(sup),
      extent_mb      = round(union_len(pmax(sup$seg_start, WIN$lo[w]),
                                       pmin(sup$seg_end,   WIN$hi[w])) / 1e6, 1),
      states         = paste(sort(unique(sup$state)), collapse = ","),
      max_llr        = if (nrow(sup)) round(max(sup$LLR, na.rm = TRUE), 1) else NA_real_,
      max_phi        = if (nrow(sup)) round(max(sup$phi_mle, na.rm = TRUE), 3) else NA_real_,
      n_genes        = if (nrow(sup)) sum(sup$n_genes, na.rm = TRUE) else NA_integer_,
      cell_frac      = round(cellfrac, 3),
      rounds_calling = ncall,
      rounds_total   = nr
    )
  }
}

res <- rbindlist(rows)

# Confidence tier. Among the 29 called rows the segment LLRs fall into two
# groups with nothing between them: 5.1-9.1 (n = 6) and >= 52 (n = 23). The
# per-cell posterior then splits the low-LLR group -- four of the six still
# carry the gain in 59-97% of cells, so the weak segment statistic is not the
# whole story; two do not.
#
#   strong    max_llr >= 50 and cell_frac >= 0.5
#   moderate  exactly one of those holds
#   weak      neither: low segment evidence AND a minority of cells
res[, support := fifelse(!called, NA_character_,
      fifelse(max_llr >= 50 & cell_frac >= 0.5, "strong",
      fifelse(max_llr >= 50 | cell_frac >= 0.5, "moderate", "weak")))]

setorder(res, event, -called, -arm_frac)
fwrite(res, "results/rb_2p_6p_gain_by_sample.csv")
fwrite(res[called == TRUE], "results/rb_2p_6p_gain_legitimate.csv")

cat("\n== legitimate calls ==\n")
print(res[called == TRUE, .(sample_id, event, selected_round, arm_frac, extent_mb,
                            max_llr, cell_frac, rounds_calling, rounds_total, support)])
cat("\n== support tiers ==\n")
print(dcast(res[called == TRUE], support ~ event, fun.aggregate = length, value.var = "sample_id"))
cat("\n== not called, but some arm coverage present (near misses) ==\n")
print(res[called == FALSE & arm_frac > 0,
          .(sample_id, event, arm_frac, extent_mb, max_llr, rounds_calling, rounds_total)])
cat("\n== counts ==\n")
print(res[, .(n_called = sum(called), n_samples = .N), by = event])
cat("\nboth 2p and 6p:",
    paste(sort(res[called == TRUE, .N, by = sample_id][N == 2, sample_id]), collapse = ", "), "\n")

# The pipeline still routes the 2p/6p collages and integrations off a
# hand-curated list in R/pipeline_targets_inputs.R (rb_scna_samples), written
# before the round revamp. Report the drift rather than silently disagreeing.
HARD <- list(
  "2p_gain" = c("SRX10264523", "SRX10264524", "SRX10264525", "SRX10264526",
                "SRX14116947", "SRX14116944"),
  "6p_gain" = c("SRX10264524", "SRX10264525", "SRX14116944")
)
cat("\n== vs hardcoded rb_scna_samples (R/pipeline_targets_inputs.R:134-135) ==\n")
for (ev in names(HARD)) {
  now <- sort(res[called == TRUE & event == ev, sample_id])
  cat("\n", ev, "\n", sep = "")
  cat("  hardcoded  n=", length(HARD[[ev]]), ": ", paste(sort(HARD[[ev]]), collapse = " "), "\n", sep = "")
  cat("  majority   n=", length(now), ": ", paste(now, collapse = " "), "\n", sep = "")
  cat("  added by revamp : ", paste(setdiff(now, HARD[[ev]]), collapse = " "), "\n", sep = "")
  cat("  dropped         : ", paste(setdiff(HARD[[ev]], now), collapse = " "), "\n", sep = "")
}
