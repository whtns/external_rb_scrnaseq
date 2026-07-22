# Tier 1: recover 6p gains numbat did not call, SRX samples only. Read-only.
#
# CORRECTED. A first version of this script ranked on `cnv_state_post` from
# bulk_clones_final.tsv.gz, treating "cnv_state_post == neu" as "not called".
# That is wrong and produced almost entirely false positives -- e.g. SRX10031193
# seg 6b showed cnv_state_post = neu with LLR 1730.7, while segs_consensus_4.tsv
# calls that very segment `amp` at LLR 561.4 and the `members` column assigns it
# to clones 3-8. bulk_clones_final's per-clone state fields are a pseudobulk
# retest, NOT the final call.
#
# Authoritative sources, used here:
#   - final segs_consensus_<N>.tsv  -> which segments were called, and as what
#   - bulk_clones_final.tsv.gz      -> per-clone evidence (p_amp, LLR, phi_mle)
#
# So a missed gain is: a chr6p region carrying amp evidence that is NOT covered
# by any amp segment in the final consensus. That is what this script looks for.
#
# Also reports called amps spanning/abutting numbat's hardcoded MHC blacklist
# (chr6:28,510,120-33,480,577, see doc/hla_masking_6p.md), whose breakpoints
# inside that 5 Mb are undetermined.
#
# Output: results/6p_gain_candidates_srx.csv + a printed summary.
# Triage list, NOT calls -- review against the heatmaps.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tidyr)
})

NUMBAT_DIR <- "output/numbat_sridhar"
P_ARM_END  <- 58.5e6
BL_START   <- 28510120
BL_END     <- 33480577
ADJ_TOL    <- 1e6
MIN_GENES  <- 10      # a segment with a handful of genes is noise, not a lead

sample_dirs <- sort(list.dirs(NUMBAT_DIR, recursive = FALSE))
sample_dirs <- sample_dirs[grepl("/SRX", sample_dirs)]
cat("SRX sample dirs:", length(sample_dirs), "\n\n")

final_consensus <- function(d) {
  fs <- list.files(d, pattern = "^segs_consensus_[0-9]+\\.tsv$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  n <- as.integer(str_match(basename(fs), "segs_consensus_([0-9]+)\\.tsv")[, 2])
  fs[which.max(n)]
}

one_sample <- function(d) {
  sid <- basename(d)
  bf <- file.path(d, "bulk_clones_final.tsv.gz")
  cf <- final_consensus(d)
  if (!file.exists(bf) || is.null(cf)) { message("  skip ", sid, ": missing inputs"); return(NULL) }

  cons <- tryCatch(readr::read_tsv(cf, show_col_types = FALSE, progress = FALSE,
                                   col_types = readr::cols(.default = readr::col_guess(),
                                                           CHROM = readr::col_character())),
                   error = function(e) NULL)
  if (is.null(cons)) { message("  !! ", sid, ": consensus unreadable"); return(NULL) }

  # Called amp intervals on chr6p -- the ground truth for "already found".
  amp_iv <- cons %>%
    filter(CHROM == "6", cnv_state == "amp", seg_start < P_ARM_END) %>%
    select(seg_start, seg_end)

  b <- tryCatch(readr::read_tsv(bf, show_col_types = FALSE, progress = FALSE,
                                col_types = readr::cols(.default = readr::col_guess(),
                                                        CHROM = readr::col_character(),
                                                        seg = readr::col_character(),
                                                        sample = readr::col_character())),
                error = function(e) NULL)
  if (is.null(b)) { message("  !! ", sid, ": bulk unreadable"); return(NULL) }

  b <- b %>% filter(CHROM == "6", !is.na(seg), !is.na(sample), seg_start < P_ARM_END)
  if (nrow(b) == 0) return(NULL)

  # Assert, do not assume, that (clone, seg) determines the segment-level fields.
  # The previous version asserted this in a comment without checking it.
  #
  # p_amp is deliberately excluded from the key: within one segment it appears
  # both as NA and as 0 on different rows (same seg_start, cnv_state, phi_mle and
  # LLR throughout), which is a missing-vs-zero encoding quirk rather than genuine
  # ambiguity. Including it made the check fire on samples that were in fact fine.
  # It is resolved below with max(na.rm), so NA/0 collapses to 0 = no amp evidence.
  chk <- b %>% distinct(sample, seg, seg_start, seg_end, cnv_state, phi_mle, LLR, n_genes) %>%
    count(sample, seg) %>% filter(n > 1)
  if (nrow(chk) > 0) {
    message("  !! ", sid, ": (clone,seg) NOT unique for decision fields -- ",
            nrow(chk), " groups; skipping sample")
    return(NULL)
  }

  safe_max <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else max(x) }

  segs <- b %>%
    group_by(sample, seg) %>%
    summarise(seg_start = first(seg_start), seg_end = first(seg_end),
              cnv_state = first(cnv_state), phi_mle = first(phi_mle),
              LLR = first(LLR), p_amp = safe_max(p_amp),
              n_genes = first(n_genes), n_snps = first(n_snps),
              .groups = "drop") %>%
    mutate(sample_id = sid)

  # Covered by an existing amp call? Any overlap counts -- we are asking whether
  # this evidence is already represented, not whether boundaries match.
  covered <- vapply(seq_len(nrow(segs)), function(i) {
    if (nrow(amp_iv) == 0) return(FALSE)
    any(amp_iv$seg_end > segs$seg_start[i] & amp_iv$seg_start < segs$seg_end[i])
  }, logical(1))

  segs %>%
    mutate(covered_by_amp_call = covered,
           n_amp_calls_6p = nrow(amp_iv),
           spans_bl = seg_start <= BL_START & seg_end >= BL_END,
           abuts_bl = abs(seg_end - BL_START) < ADJ_TOL | abs(seg_start - BL_END) < ADJ_TOL)
}

p6 <- purrr::map(sample_dirs, one_sample) %>% purrr::compact() %>% bind_rows()
if (nrow(p6) == 0) { cat("no chr6p segments read\n"); quit(save = "no") }

p6 <- p6 %>% mutate(Mb = sprintf("%.1f-%.1f", seg_start/1e6, seg_end/1e6))
readr::write_csv(p6, "results/6p_gain_candidates_srx.csv")
cat("wrote results/6p_gain_candidates_srx.csv --", nrow(p6), "chr6p (clone,segment) records across",
    n_distinct(p6$sample_id), "SRX samples\n\n")

show <- function(df, title) {
  cat("===", title, "===\n")
  if (nrow(df) == 0) { cat("  (none)\n\n"); return(invisible()) }
  df %>%
    transmute(sample_id, clone = sample, seg, Mb,
              p_amp = round(p_amp, 4), phi = round(phi_mle, 2),
              LLR = round(LLR, 1), genes = n_genes, snps = n_snps,
              state = cnv_state, amp6p = n_amp_calls_6p,
              bl = ifelse(spans_bl, "spans", ifelse(abuts_bl, "abuts", ""))) %>%
    as.data.frame() %>% print(row.names = FALSE)
  cat("\n")
}

uncalled <- p6 %>% filter(!covered_by_amp_call, !is.na(p_amp), n_genes >= MIN_GENES)

show(uncalled %>% filter(p_amp > 0.9) %>% arrange(desc(p_amp), desc(LLR)) %>% head(40),
     paste0("A. STRONG: amp evidence with NO overlapping amp call (p_amp > 0.9, >=", MIN_GENES, " genes)"))

show(uncalled %>% filter(p_amp > 0.5, p_amp <= 0.9) %>% arrange(desc(p_amp)) %>% head(40),
     "B. SUGGESTIVE: 0.5 < p_amp <= 0.9, no overlapping amp call")

show(p6 %>% filter(covered_by_amp_call, cnv_state == "amp", spans_bl | abuts_bl) %>%
       arrange(sample_id, seg) %>% head(40),
     "C. CALLED amps spanning/abutting the MHC blacklist (breakpoints undetermined)")

cat("=== per-sample summary ===\n")
p6 %>%
  group_by(sample_id) %>%
  summarise(clones = n_distinct(sample), segs_6p = n(),
            amp_calls_6p = first(n_amp_calls_6p),
            strong = sum(!covered_by_amp_call & !is.na(p_amp) & p_amp > 0.9 & n_genes >= MIN_GENES),
            suggestive = sum(!covered_by_amp_call & !is.na(p_amp) & p_amp > 0.5 & p_amp <= 0.9 & n_genes >= MIN_GENES),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nNOTE: 'called' means an overlapping amp segment in the FINAL segs_consensus,\n")
cat("which is authoritative. bulk_clones_final's per-clone cnv_state/cnv_state_post\n")
cat("is a pseudobulk retest and must NOT be read as the call.\n")
cat("\nDONE\n")
