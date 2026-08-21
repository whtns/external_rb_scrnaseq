# Cross-round recovery of canonical RB SCNAs.
#
# numbat's consensus iteration does not converge in this cohort (only 4/71 runs
# reach d_k = 0; see src/diag_consensus_convergence.R and
# doc/diploid_baseline_pin_6p.md), so the final segs_consensus_<K>.tsv is a
# stopping point, not an answer. A segment called in rounds 1-3 and absent at
# round 4 is not evidence against the call -- SRX10264524 calls 1 canonical event
# at round 4 and all 5 across rounds 1-4.
#
# This reads every segs_consensus_<k>.tsv already on disk and scores each segment
# by how many rounds called it. No rerun.
#
# Outputs:
#   results/rb_scna_segments_by_round.csv    every non-neutral segment x sample,
#                                            with n_rounds_called / n_rounds,
#                                            max LLR, and called_final
#   results/rb_scna_canonical_by_sample.csv  sample x canonical event, final vs
#                                            union, with stability and max LLR
#
# Usage:
#   Rscript src/diag_rb_scna_recovery.R [numbat_dir] [suffix]
#   e.g.    Rscript src/diag_rb_scna_recovery.R output/numbat_multiallelic _ma
#
# The five canonical RB arm events. hg38 centromere midpoints; chr13 is
# acrocentric so all of it is 13q for this purpose.
#   1q  gain   CHROM 1,  seg_end   > 125.0 Mb
#   2p  gain   CHROM 2,  seg_start <  93.0 Mb   (MYCN at 15.9 Mb)
#   6p  gain   CHROM 6,  seg_start <  59.0 Mb
#   13q loss   CHROM 13                          (RB1 at 48.3 Mb)
#   16q loss   CHROM 16, seg_end   >  36.8 Mb

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tidyr)
})

args       <- commandArgs(trailingOnly = TRUE)
NUMBAT_DIR <- if (length(args) >= 1) args[[1]] else "output/numbat_sridhar"
SUFFIX     <- if (length(args) >= 2) args[[2]] else ""

GAIN <- "amp"        # matches amp and bamp
LOSS <- "del|loh"    # matches del, bdel, loh

classify_event <- function(CHROM, seg_start, seg_end, state) {
  gain <- str_detect(state, GAIN)
  loss <- str_detect(state, LOSS)
  case_when(
    CHROM == "1"  & seg_end    > 125e6  & gain ~ "1q_gain",
    CHROM == "2"  & seg_start  <  93e6  & gain ~ "2p_gain",
    CHROM == "6"  & seg_start  <  59e6  & gain ~ "6p_gain",
    CHROM == "13"                       & loss ~ "13q_loss",
    CHROM == "16" & seg_end    > 36.8e6 & loss ~ "16q_loss",
    TRUE ~ NA_character_
  )
}

read_round <- function(f) {
  x <- tryCatch(
    read_tsv(f, show_col_types = FALSE, progress = FALSE,
             col_types = cols(.default = col_guess(), CHROM = col_character())),
    error = function(e) NULL
  )
  if (is.null(x) || !all(c("CHROM", "seg_start", "seg_end", "cnv_state") %in% names(x))) {
    return(NULL)
  }
  # cnv_state_post is the multi-allelic-aware call where present; fall back to
  # cnv_state, which is all the multi_allelic = FALSE runs ever wrote.
  state <- if ("cnv_state_post" %in% names(x)) {
    coalesce(na_if(as.character(x$cnv_state_post), ""), as.character(x$cnv_state))
  } else {
    as.character(x$cnv_state)
  }
  if (!"LLR" %in% names(x)) x$LLR <- NA_real_
  x %>%
    mutate(
      state = state,
      round = as.integer(str_match(basename(f), "segs_consensus_([0-9]+)\\.tsv")[, 2]),
      # File mtime, to detect dirs whose "rounds" are spliced from different runs.
      mtime = file.info(f)$mtime
    ) %>%
    filter(!is.na(state), state != "neu") %>%
    select(round, CHROM, seg, seg_start, seg_end, state, LLR, mtime)
}

sample_dirs <- list.dirs(NUMBAT_DIR, recursive = FALSE)
message("scanning ", length(sample_dirs), " sample dirs under ", NUMBAT_DIR)

per_sample <- map(sample_dirs, function(d) {
  fs <- list.files(d, pattern = "^segs_consensus_[0-9]+\\.tsv$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  rounds <- map_dfr(fs, read_round)
  if (nrow(rounds) == 0) return(NULL)
  rounds %>% mutate(sample_id = basename(d), n_rounds = length(fs))
})

segs <- bind_rows(per_sample)
stopifnot(nrow(segs) > 0)

# PROVENANCE GUARD.
# The cross-round union assumes segs_consensus_1..K in a dir are successive
# rounds of ONE run. That is not always true: some dirs hold rounds spliced from
# different runs (a later rerun crashed partway and overwrote only the early
# rounds, leaving older files behind). In those dirs a round-to-round difference
# is partly a run-to-run difference under possibly different parameters, so the
# stability score is not interpretable. Detected by mtime spread across rounds.
MAX_MTIME_SPAN_DAYS <- 1

provenance <- segs %>%
  group_by(sample_id) %>%
  summarise(
    mtime_span_days = as.numeric(difftime(max(mtime), min(mtime), units = "days")),
    .groups = "drop"
  ) %>%
  mutate(mixed_provenance = mtime_span_days > MAX_MTIME_SPAN_DAYS)

segs <- segs %>% left_join(provenance, by = "sample_id")

if (any(provenance$mixed_provenance)) {
  message("mixed-provenance sample dirs (rounds not from one run): ",
          paste(provenance$sample_id[provenance$mixed_provenance], collapse = ", "))
}

max_or_na <- function(x) {
  v <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.finite(v)) v else NA_real_
}

# One row per (sample, segment interval, state), scored across rounds. Keyed on
# the interval and state rather than numbat's `seg` label, which is re-lettered
# every round and so is not comparable between them.
segments <- segs %>%
  group_by(sample_id, n_rounds, mixed_provenance, CHROM, seg_start, seg_end, state) %>%
  summarise(
    n_rounds_called = n_distinct(round),
    rounds_called   = paste(sort(unique(round)), collapse = ","),
    max_LLR         = max_or_na(LLR),
    called_final    = any(round == first(n_rounds)),
    .groups = "drop"
  ) %>%
  mutate(
    stability = n_rounds_called / n_rounds,
    event     = classify_event(CHROM, seg_start, seg_end, state)
  ) %>%
  arrange(sample_id, CHROM, seg_start)

# Per sample x canonical event: is it called at the final round, and is it called
# anywhere at all?
canonical <- segments %>%
  filter(!is.na(event)) %>%
  group_by(sample_id, event, mixed_provenance) %>%
  summarise(
    called_final   = any(called_final),
    best_stability = max(stability),
    max_LLR        = max_or_na(max_LLR),
    rounds_called  = paste(sort(unique(as.integer(unlist(str_split(rounds_called, ","))))),
                           collapse = ","),
    n_rounds       = first(n_rounds),
    .groups = "drop"
  ) %>%
  mutate(
    recovered_by_union = !called_final,
    # A call seen in a single round of a 4-round run is weaker evidence than one
    # seen in three. Tier accordingly rather than treating the union as flat.
    tier = case_when(
      called_final          ~ "called",
      best_stability >= 0.5 ~ "recovered_stable",
      TRUE                  ~ "recovered_marginal"
    )
  ) %>%
  arrange(sample_id, event)

dir.create("results", showWarnings = FALSE)
seg_out <- file.path("results", paste0("rb_scna_segments_by_round", SUFFIX, ".csv"))
can_out <- file.path("results", paste0("rb_scna_canonical_by_sample", SUFFIX, ".csv"))
write_csv(segments,  seg_out)
write_csv(canonical, can_out)

n_final <- sum(canonical$called_final)
n_union <- nrow(canonical)

clean   <- canonical %>% filter(!mixed_provenance)
c_final <- sum(clean$called_final)
c_union <- nrow(clean)

cat("\n=== canonical RB events,", basename(NUMBAT_DIR), "===\n")
cat("samples:", n_distinct(segments$sample_id), "\n")
cat(sprintf("ALL dirs     final: %d   union: %d   (+%d, %.0f%%)\n",
            n_final, n_union, n_union - n_final,
            100 * (n_union - n_final) / n_final))
cat(sprintf("CLEAN dirs   final: %d   union: %d   (+%d, %.0f%%)   <- the number to quote\n",
            c_final, c_union, c_union - c_final,
            100 * (c_union - c_final) / c_final))
cat(sprintf("excluded as mixed-provenance: %d samples, %d events\n",
            n_distinct(canonical$sample_id[canonical$mixed_provenance]),
            sum(canonical$mixed_provenance)))

cat("\nby event:\n")
canonical %>%
  group_by(event) %>%
  summarise(final = sum(called_final), union = n(),
            recovered_stable   = sum(tier == "recovered_stable"),
            recovered_marginal = sum(tier == "recovered_marginal"),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nsamples gaining >=1 event from the union:",
    n_distinct(canonical$sample_id[!canonical$called_final]),
    "of", n_distinct(canonical$sample_id), "\n")

cat("\nrecovered events by tier (clean dirs only):\n")
clean %>% filter(!called_final) %>% count(tier) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nwrote", seg_out, "and", can_out, "\n")
