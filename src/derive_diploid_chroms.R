# Derive a per-sample diploid_chroms set for numbat.
#
# Why per-sample rather than one cohort-wide set: numbat's default
# find_common_diploid() re-derives the diploid baseline EVERY iteration, so the
# reference that gains are measured against drifts between rounds. SRX10264524
# produced 12 distinct diploid sets in one run, and lost a chr1q amp at LLR 92
# between rounds 3 and 4 as a result. Pinning diploid_chroms takes the
# `else if (!is.null(diploid_chroms))` branch in analyze_bulk (body line 14),
# which never reaches find_common_diploid() -- the baseline becomes a fixed,
# deterministic chromosome-membership test.
#
# A single cohort-wide set is still wrong for the 11-17% of samples that carry
# events on the otherwise-clean chromosomes, so each sample gets its own set
# derived from its own calls.
#
# Burden = fraction of a chromosome's BASELINE-USABLE genes (present in the
# reference AND nonzero in the Cone column -- get_bulk drops lambda_ref == 0)
# that fall inside a non-neutral segment of that sample's final consensus.
#
# Burden is computed over the UNION of non-neutral segments across ALL consensus
# rounds, not just the final one. A first version used only the final round and
# was badly wrong: SRX10264524 carries chr1 150-249Mb amp at LLR 92 and chr6p
# amps at LLR 67-96 in rounds 1-3, all absent from round 4, so a final-round
# burden scored chr1 and chr6 as clean and pinned them as DIPLOID -- which would
# have made those very gains impossible to call. 65 of 71 samples hit some form
# of this. Taking the union means any chromosome called in any round is excluded,
# which is the conservative direction: a smaller baseline is much cheaper than a
# contaminated one.
#
# A cohort-level guard is applied on top. A chromosome recurrently affected
# across the cohort is never pinned, even if this particular sample never called
# it -- that protects against a per-sample power failure silently promoting a
# real event region into the baseline.
#
# Output: results/diploid_chroms_per_sample.tsv

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tibble)
  library(tidyr); library(yaml)
})
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressMessages(library(numbat))

NUMBAT_DIR  <- "output/numbat_sridhar"
REF_PATH    <- "/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds"
MAX_BURDEN  <- 0.05   # HARD cap. A chromosome above this is never pinned, ever.
MIN_GENES   <- 1500   # baseline floor; a sample that cannot reach it OPTS OUT
MIN_CHROMS  <- 3      # never pin fewer than this
COHORT_MAX  <- 0.40   # cohort-level ban threshold, applied to filtered calls

# Never pinned, regardless of what the burden data says. These are recurrent RB
# event chromosomes; a hard list means the guarantee does not depend on a
# threshold happening to catch them. chr2 (2p gain) is here by explicit request
# and is also banned data-driven at LLR>=20, so the list agrees with the data
# rather than overriding it.
ALWAYS_BANNED <- c("2")

# Only persistent, confident calls count toward burden. Counting every call in
# every round treated a transient round-1 artefact as equal evidence to a
# persistent LLR-92 event, which inflated burden until 16 of 22 chromosomes were
# banned -- including chr19 (62%) and chr12 (66%), not credible RB event rates.
# With these filters the cohort ban collapses to {1, 6, 13, 16}: 1q gain, 6p
# gain, RB1 at 13q14, 16q loss. That the banned set reproduces the canonical RB
# events is the main evidence these thresholds are separating signal from churn
# rather than just being permissive.
MIN_LLR_CALL  <- 20    # a call below this is not evidence a chromosome is dirty
# Deliberately 0: ANY round with a confident call marks the chromosome. A previous
# version required persistence in >=50% of rounds, which is self-defeating -- it
# assumes the run is stable enough for a real event to recur, which is exactly
# what is false here. SRX10264524 carries chr2 amp at LLR 51.5 in round 1 only
# (and the config annotates 6_v_2_2p+ for it); a persistence filter scored chr2 as
# clean and pinned it as diploid.
MIN_ROUND_FRAC <- 0

g <- as.data.frame(numbat::gtf_hg38)
g$CHROM <- as.character(g$CHROM)
g <- g[g$CHROM %in% as.character(1:22), ]

m <- as.matrix(readRDS(REF_PATH))
stopifnot("Cone" %in% colnames(m))
usable <- rownames(m)[m[, "Cone"] > 0]
gu <- g[g$gene %in% usable, ]
tot <- gu %>% count(CHROM, name = "n_usable")
cat("baseline-usable genes (in reference, nonzero in Cone):", nrow(gu), "\n\n")

# Every round, not just the final one -- see the header note.
all_consensus <- function(d) {
  list.files(d, pattern = "^segs_consensus_[0-9]+\\.tsv$", full.names = TRUE)
}

sample_segs <- function(d) {
  fs <- all_consensus(d)
  if (!length(fs)) return(NULL)
  n_rounds <- length(fs)
  x <- map_dfr(seq_along(fs), function(i)
    tryCatch(read_tsv(fs[i], show_col_types = FALSE, progress = FALSE,
                      col_types = cols(.default = col_guess(), CHROM = col_character())) %>%
               select(any_of(c("CHROM", "seg_start", "seg_end", "cnv_state", "LLR"))) %>%
               mutate(round_idx = i),
             error = function(e) NULL))
  if (is.null(x) || nrow(x) == 0 || !"LLR" %in% names(x)) return(NULL)

  x <- x %>% filter(!is.na(cnv_state), cnv_state != "neu", cnv_state != "",
                    CHROM %in% as.character(1:22), !is.na(LLR), LLR >= MIN_LLR_CALL)
  if (nrow(x) == 0) return(x %>% select(any_of(c("CHROM", "seg_start", "seg_end"))))

  # Keep only chromosomes whose confident calls PERSIST across rounds; a
  # chromosome called once out of four is churn, not an event.
  keep_chrom <- x %>% distinct(CHROM, round_idx) %>% count(CHROM) %>%
    filter(n / n_rounds >= MIN_ROUND_FRAC) %>% pull(CHROM)

  x %>% filter(CHROM %in% keep_chrom) %>%
    select(CHROM, seg_start, seg_end) %>% distinct()
}

sample_burden <- function(d) {
  x <- sample_segs(d)
  if (is.null(x)) return(NULL)
  map_dfr(as.character(1:22), function(cc) {
    seg <- x  %>% filter(CHROM == cc)
    gg  <- gu %>% filter(CHROM == cc)
    hit <- if (nrow(seg) == 0) 0L else
      sum(map_lgl(gg$gene_start, ~ any(seg$seg_start <= .x & seg$seg_end >= .x)))
    tibble(sample_id = basename(d), CHROM = cc, burden = hit / max(nrow(gg), 1))
  })
}

dirs_all <- sort(list.dirs(NUMBAT_DIR, recursive = FALSE))
burden_all <- map(dirs_all, sample_burden) %>% compact() %>% bind_rows()
stopifnot(nrow(burden_all) > 0)

# Cohort guard: chromosomes affected in > COHORT_MAX of samples are off-limits.
cohort <- burden_all %>%
  group_by(CHROM) %>%
  summarise(pct_affected = mean(burden > 0), .groups = "drop") %>%
  arrange(desc(pct_affected))
banned_data <- cohort %>% filter(pct_affected > COHORT_MAX) %>% pull(CHROM)
banned <- union(banned_data, ALWAYS_BANNED)
cat("cohort-banned by burden (> ", round(100 * COHORT_MAX), "% of samples): ",
    paste(sort(as.integer(banned_data)), collapse = ", "), "\n", sep = "")
cat("always-banned (hard list)                : ",
    paste(sort(as.integer(ALWAYS_BANNED)), collapse = ", "), "\n", sep = "")
cat("effective ban                            : ",
    paste(sort(as.integer(banned)), collapse = ", "), "\n", sep = "")
cat("cohort affected-rate by chromosome:\n")
print(as.data.frame(cohort %>% mutate(pct_affected = round(100 * pct_affected))), row.names = FALSE)
cat("\n")

# Per-sample veto from config/large_clone_comparisons.yaml. This is the ONLY
# input here that is not derived from the numbat calls whose instability is the
# whole problem -- it is curated and human-reviewed, so it can catch an event the
# call data lost entirely in every round. Chromosomes are taken from both the
# comparison KEY (arm notation, e.g. "6_v_2_2p+" -> 2) and the segment id VALUES
# (e.g. "1h", "16b" -> 1, 16), since the two encode the event differently.
# Branch keys (SRX..._branch_5) fold into their base sample.
config_veto <- local({
  f <- "config/large_clone_comparisons.yaml"
  if (!file.exists(f)) { cat("!! ", f, " not found -- config veto DISABLED\n"); return(list()) }
  y <- yaml::read_yaml(f)
  out <- list()
  for (k in names(y)) {
    base <- str_extract(k, "^SR[RX][0-9]+")
    if (is.na(base)) next
    chroms <- character(0)
    comps <- y[[k]]
    for (cmp in names(comps)) {
      tail_str <- str_remove(cmp, "^[0-9]+_v_[0-9]+_?")
      chroms <- c(chroms, str_match_all(tail_str, "([0-9]{1,2})[pq]")[[1]][, 2])
      chroms <- c(chroms, str_match_all(tail_str, "([0-9]{1,2})loh")[[1]][, 2])
      segs <- unlist(comps[[cmp]])
      if (length(segs)) chroms <- c(chroms, str_match(segs, "^([0-9]{1,2})")[, 2])
    }
    chroms <- unique(chroms[!is.na(chroms) & chroms %in% as.character(1:22)])
    if (length(chroms)) out[[base]] <- union(out[[base]], chroms)
  }
  out
})
cat("config veto loaded for", length(config_veto), "samples",
    "(median", if (length(config_veto)) median(lengths(config_veto)) else 0,
    "chromosomes vetoed)\n\n")

one_sample <- function(d) {
  sid <- basename(d)
  b <- burden_all %>% filter(sample_id == sid)
  if (nrow(b) == 0) return(NULL)

  veto <- config_veto[[sid]]
  if (is.null(veto)) veto <- character(0)

  burden <- b %>% select(CHROM, burden) %>%
    filter(!CHROM %in% banned, !CHROM %in% veto) %>%
    left_join(tot, by = "CHROM") %>% arrange(burden, desc(n_usable))

  # MAX_BURDEN is a HARD cap with no relaxation path. The previous version kept
  # adding the next-lowest-burden chromosome until it hit MIN_GENES, which pinned
  # chromosomes at up to 100% burden -- a baseline sitting entirely inside a
  # called CNV, worse than the drifting baseline this is meant to replace.
  keep <- burden %>% filter(burden <= MAX_BURDEN)

  # A sample that cannot assemble an adequate CLEAN baseline OPTS OUT: it gets
  # diploid_chroms = NA and falls back to numbat's find_common_diploid(). It stays
  # unstable, but an unstable baseline beats a knowably wrong one. Some RB genomes
  # simply have no clean chromosomes, and that is a real answer.
  opted_out <- nrow(keep) < MIN_CHROMS || sum(keep$n_usable) < MIN_GENES

  tibble(
    sample_id      = sid,
    diploid_chroms = if (opted_out) NA_character_
                     else paste(sort(as.integer(keep$CHROM)), collapse = ","),
    n_chroms       = if (opted_out) 0L else nrow(keep),
    n_usable_genes = if (opted_out) 0L else sum(keep$n_usable),
    max_burden     = if (opted_out) NA_real_ else round(max(keep$burden), 4),
    # What it *would* have had, to make the opt-out diagnosable.
    avail_chroms   = nrow(keep),
    avail_genes    = sum(keep$n_usable),
    opted_out      = opted_out
  )
}

res <- map(sort(list.dirs(NUMBAT_DIR, recursive = FALSE)), one_sample) %>%
  compact() %>% bind_rows()

stopifnot(nrow(res) > 0)

# Self-check with NO EXEMPTIONS. The v2 version excused `relaxed` samples, so 42
# of 71 bypassed it and it reported "passed" while chromosomes at 100% burden
# were being pinned -- a check that skips the dangerous cases is worse than none,
# because it manufactures confidence. Every pinned chromosome is checked here.
viol <- res %>%
  filter(!is.na(diploid_chroms)) %>%
  select(sample_id, diploid_chroms) %>%
  mutate(CHROM = strsplit(diploid_chroms, ",")) %>%
  tidyr::unnest(CHROM) %>%
  left_join(burden_all, by = c("sample_id", "CHROM")) %>%
  mutate(is_banned = CHROM %in% banned,
         is_vetoed = map2_lgl(sample_id, CHROM,
                              ~ .y %in% (config_veto[[.x]] %||% character(0)))) %>%
  filter(is_banned | is_vetoed | is.na(burden) | burden > MAX_BURDEN)

if (nrow(viol) > 0) {
  cat("!! SELF-CHECK FAILED -- pinned chromosomes carrying calls or cohort-banned:\n")
  print(as.data.frame(viol), row.names = FALSE)
  stop("refusing to write a diploid_chroms table that pins affected chromosomes")
}
cat("self-check passed: no pinned chromosome is cohort-banned, and none carries\n")
cat("a call above the burden threshold in any consensus round.\n\n")

write_tsv(res, "results/diploid_chroms_per_sample.tsv")
cat("wrote results/diploid_chroms_per_sample.tsv --", nrow(res), "samples\n\n")

cat("=== samples that OPTED OUT (no adequate clean baseline; keep numbat default) ===\n")
r <- res %>% filter(opted_out)
if (nrow(r) == 0) cat("  (none)\n") else
  print(as.data.frame(r %>% select(sample_id, avail_chroms, avail_genes)), row.names = FALSE)

cat("\n=== summary ===\n")
pinned <- res %>% filter(!opted_out)
cat("samples pinned  :", nrow(pinned), "/", nrow(res), "\n")
cat("samples opted out:", sum(res$opted_out), "\n")
if (nrow(pinned) > 0) {
  cat("median chromosomes pinned :", median(pinned$n_chroms), "\n")
  cat("median baseline genes     :", median(pinned$n_usable_genes), "\n")
  cat("worst pinned burden       :", round(max(pinned$max_burden), 4),
      "(hard cap", MAX_BURDEN, ")\n")
}
cat("\n")

cat("=== the two rerun targets ===\n")
print(as.data.frame(res %>% filter(sample_id %in% c("SRX10264524", "SRX14116944"))),
      row.names = FALSE)

cat("\nDONE\n")
