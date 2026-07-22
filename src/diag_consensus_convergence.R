# How many numbat iterations were actually needed?
#
# Answers it empirically rather than by rule of thumb, using data already on
# disk: numbat writes segs_consensus_<k>.tsv for every iteration, so the change
# between consecutive rounds is recoverable for every sample with no rerun.
#
# Metric: d_k = size of the symmetric difference between the NON-NEUTRAL segment
# sets at round k and round k-1, keyed on (CHROM, seg_start, seg_end, cnv_state).
# d_k = 0 means that iteration changed nothing, i.e. the run had settled.
#
# The decision this supports:
#   max_iter = smallest K with d_K = 0 for the large majority of samples
#   samples still changing at K are OSCILLATING, not under-iterated -- marginal
#   calls flickering across min_LLR. More rounds alternate the answer rather than
#   improving it, so they need a stability flag or a wider min_LLR margin.
#
# `oscillating` below is detected directly: round k identical to round k-2 but
# different from k-1 is a 2-cycle, which no amount of extra iterations resolves.
#
# Runs over ALL samples, not just SRX -- the files are small TSVs.
# Output: results/consensus_convergence.csv + a printed cohort summary.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tidyr)
})

NUMBAT_DIR <- "output/numbat_sridhar"

seg_key <- function(f) {
  x <- tryCatch(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE,
                                col_types = readr::cols(.default = readr::col_guess(),
                                                        CHROM = readr::col_character())),
                error = function(e) NULL)
  if (is.null(x) || !all(c("CHROM", "seg_start", "seg_end", "cnv_state") %in% names(x))) return(NULL)
  x %>%
    filter(!is.na(cnv_state), cnv_state != "neu", cnv_state != "") %>%
    transmute(k = paste(CHROM, seg_start, seg_end, cnv_state, sep = ":")) %>%
    pull(k) %>% unique() %>% sort()
}

one_sample <- function(d) {
  sid <- basename(d)
  fs <- list.files(d, pattern = "^segs_consensus_[0-9]+\\.tsv$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  n <- as.integer(str_match(basename(fs), "segs_consensus_([0-9]+)\\.tsv")[, 2])
  fs <- fs[order(n)]; n <- sort(n)

  sets <- purrr::map(fs, seg_key)
  ok <- !purrr::map_lgl(sets, is.null)
  sets <- sets[ok]; n <- n[ok]
  if (length(sets) == 0) return(NULL)

  purrr::map_dfr(seq_along(sets), function(i) {
    prev  <- if (i > 1) sets[[i - 1]] else NULL
    prev2 <- if (i > 2) sets[[i - 2]] else NULL
    tibble::tibble(
      sample_id   = sid,
      round       = n[i],
      n_segs      = length(sets[[i]]),
      d_prev      = if (is.null(prev)) NA_integer_
                    else length(union(setdiff(sets[[i]], prev), setdiff(prev, sets[[i]]))),
      # A 2-cycle: identical to two rounds back but not to the previous round.
      oscillating = !is.null(prev2) &&
                    setequal(sets[[i]], prev2) && !setequal(sets[[i]], prev)
    )
  })
}

res <- purrr::map(sort(list.dirs(NUMBAT_DIR, recursive = FALSE)), one_sample) %>%
  purrr::compact() %>% bind_rows()
if (nrow(res) == 0) { cat("no consensus files found\n"); quit(save = "no") }

readr::write_csv(res, "results/consensus_convergence.csv")
cat("wrote results/consensus_convergence.csv --", n_distinct(res$sample_id), "samples\n\n")

cat("=== how many rounds each run actually got ===\n")
res %>% group_by(sample_id) %>% summarise(max_round = max(round), .groups = "drop") %>%
  count(max_round, name = "n_samples") %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n=== change per round (d_k = segments differing from previous round) ===\n")
res %>% filter(!is.na(d_prev)) %>%
  group_by(round) %>%
  summarise(n_samples = n(),
            settled = sum(d_prev == 0),
            pct_settled = round(100 * mean(d_prev == 0)),
            median_d = median(d_prev), max_d = max(d_prev), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# The decision: per sample, the first round from which nothing changes again.
per_sample <- res %>%
  group_by(sample_id) %>%
  arrange(round, .by_group = TRUE) %>%
  summarise(
    max_round   = max(round),
    last_d      = last(d_prev),
    ever_osc    = any(oscillating, na.rm = TRUE),
    settled_at  = {
      dd <- d_prev; rr <- round
      idx <- which(!is.na(dd))
      hit <- NA_integer_
      for (j in idx) if (all(dd[idx[idx >= j]] == 0)) { hit <- rr[j]; break }
      hit
    },
    .groups = "drop")

cat("\n=== per-sample outcome ===\n")
cat("settled (no further change through the last round):",
    sum(!is.na(per_sample$settled_at)), "/", nrow(per_sample), "\n")
cat("still changing at the last round               :",
    sum(is.na(per_sample$settled_at)), "\n")
cat("showing a 2-cycle (oscillating)                :", sum(per_sample$ever_osc), "\n\n")

cat("If you had stopped at round K, what fraction were settled by then?\n")
for (K in sort(unique(res$round))) {
  elig <- per_sample %>% filter(max_round >= K)
  if (nrow(elig) == 0) next
  good <- sum(!is.na(elig$settled_at) & elig$settled_at <= K)
  cat(sprintf("  K = %d : %3d / %3d settled (%3.0f%%)\n", K, good, nrow(elig),
              100 * good / nrow(elig)))
}

cat("\n=== samples STILL CHANGING at their last round (more iterations may help) ===\n")
still <- per_sample %>% filter(is.na(settled_at), !ever_osc) %>% arrange(desc(last_d))
if (nrow(still) == 0) cat("  (none)\n") else
  still %>% select(sample_id, max_round, last_d) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n=== OSCILLATING samples (more iterations will NOT help) ===\n")
osc <- per_sample %>% filter(ever_osc) %>% arrange(desc(last_d))
if (nrow(osc) == 0) cat("  (none)\n") else
  osc %>% select(sample_id, max_round, last_d) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nDONE\n")
