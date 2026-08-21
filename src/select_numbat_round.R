# Pick and record a favorable numbat consensus round per sample.
#
# WHY A ROUND, NOT THE UNION
#
# numbat's consensus iteration does not converge in this cohort, and it is not
# monotone: canonical RB events called per round run 106 / 98 / 92 / 86 across
# the 33 clean 4-round samples (src/diag_rb_scna_recovery.R). Production reports
# the LAST round, which is the worst one. The cross-round union recovers +37
# events, but a union is NOT a numbat state -- it merges segmentations from
# different iterations that disagree about segment boundaries, so it cannot
# produce a clone tree, a heatmap, or anything the downstream pipeline consumes.
# It is a scoring device, not an object.
#
# A selected round is a real numbat state. Every per-round artefact is mutually
# consistent for a given k: segs_consensus_k, clone_post_k, bulk_clones_k,
# exp_post_k, allele_post_k, joint_post_k, clones_k.rds, tree_final_k.rds,
# mut_graph_k.rds, geno_k.tsv. Selecting k therefore yields something usable
# downstream, which the union never can.
#
# WHAT "FAVORABLE" MEANS HERE -- AND THE BIAS THIS CARRIES
#
# Read this before quoting any number out of this script.
#
# Choosing the round that maximises canonical RB events is selecting on the
# outcome. The resulting per-sample count is a MAXIMUM over rounds, not an
# unbiased estimate, and it is inflated relative to any single pre-specified
# round. That is acceptable for the stated purpose -- use the RB prior to avoid
# reporting an arbitrarily bad iteration -- but the inflation is real and the
# selected-round total must never be compared against a fixed-round total as
# though both were unbiased.
#
# There is also a specific failure mode to watch. Earlier rounds call MORE of
# everything (830 non-neutral segments at round 1 vs 642 at round 4) while the
# canonical FRACTION stays pinned at ~26%. So "most canonical events" can
# degenerate into "least filtered round". Two guards:
#
#   1. The default criterion counts only CORROBORATED canonical events -- those
#      the sample calls in at least MIN_ROUNDS rounds. A singleton call that
#      appears in one round and nowhere else does not earn the selection.
#   2. n_segs is reported for the selected and final round, and the summary
#      prints the distribution of selected round indices. If selection collapses
#      onto round 1 everywhere, the honest description is "use round 1", not
#      "per-sample selection", and the output makes that visible.
#
# Both criteria are computed and written out, so the raw-max alternative can be
# inspected rather than taken on faith.
#
# Usage:
#   Rscript src/select_numbat_round.R [numbat_dir] [suffix]
#
# Outputs:
#   results/numbat_selected_round[suffix].csv    one row per sample: chosen k,
#                                                per-round scores, why
#   results/numbat_selected_paths[suffix].csv    resolved per-round file paths
#                                                for the chosen k

suppressPackageStartupMessages({
  library(data.table)
})

args       <- commandArgs(trailingOnly = TRUE)
NUMBAT_DIR <- if (length(args) >= 1) args[[1]] else "output/numbat_sridhar"
SUFFIX     <- if (length(args) >= 2) args[[2]] else ""

MIN_ROUNDS <- 2     # rounds an event must appear in to count as corroborated
MAX_MTIME_SPAN_DAYS <- 1

# The five canonical RB arm events, matching src/diag_rb_scna_recovery.R.
GAIN <- "amp"
LOSS <- "del|loh"

classify_event <- function(CHROM, seg_start, seg_end, state) {
  gain <- grepl(GAIN, state)
  loss <- grepl(LOSS, state)
  fifelse(CHROM == "1"  & seg_end   > 125e6  & gain, "1q_gain",
  fifelse(CHROM == "2"  & seg_start <  93e6  & gain, "2p_gain",
  fifelse(CHROM == "6"  & seg_start <  59e6  & gain, "6p_gain",
  fifelse(CHROM == "13"                      & loss, "13q_loss",
  fifelse(CHROM == "16" & seg_end   > 36.8e6 & loss, "16q_loss",
          NA_character_)))))
}

# Per-round artefact families. A selected round resolves to exactly these.
# A round is only usable downstream if its clone/tree artefacts were actually
# written. Two samples (SRR13633759, SRX10031191) died partway through their
# last round: segs_consensus_2.tsv exists but clones_2.rds and tree_final_2.rds
# do not. Selecting such a round yields something that cannot build a clone tree,
# so completeness is a hard constraint -- an incomplete round is chosen only if
# no complete round exists for that sample.
REQUIRED_FILES <- c(
  "segs_consensus_%d.tsv", "clone_post_%d.tsv",
  "clones_%d.rds", "tree_final_%d.rds", "mut_graph_%d.rds"
)

ROUND_FILES <- c(
  "segs_consensus_%d.tsv", "clone_post_%d.tsv", "exp_post_%d.tsv",
  "allele_post_%d.tsv", "joint_post_%d.tsv", "geno_%d.tsv",
  "bulk_clones_%d.tsv.gz", "bulk_subtrees_%d.tsv.gz",
  "clones_%d.rds", "subtrees_%d.rds", "tree_final_%d.rds",
  "treeML_%d.rds", "tree_list_%d.rds", "mut_graph_%d.rds"
)

read_rounds <- function(dir) {
  fs <- list.files(dir, pattern = "^segs_consensus_[0-9]+\\.tsv$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  rbindlist(lapply(fs, function(f) {
    z <- tryCatch(fread(f, showProgress = FALSE, colClasses = list(character = "CHROM")),
                  error = function(e) NULL)
    if (is.null(z) || nrow(z) == 0) return(NULL)
    st <- if ("cnv_state_post" %in% names(z)) z$cnv_state_post else z$cnv_state
    data.table(
      round = as.integer(sub(".*_([0-9]+)\\.tsv$", "\\1", basename(f))),
      CHROM = as.character(z$CHROM),
      seg_start = as.numeric(z$seg_start), seg_end = as.numeric(z$seg_end),
      state = st,
      LLR   = suppressWarnings(as.numeric(z$LLR)),
      mtime = file.info(f)$mtime
    )
  }), fill = TRUE)
}

dirs <- list.dirs(NUMBAT_DIR, recursive = FALSE)
dirs <- dirs[file.exists(file.path(dirs, "segs_consensus_1.tsv"))]
cat("sample dirs with consensus rounds:", length(dirs), "\n")
if (length(dirs) == 0) stop("nothing to do")

rows <- list(); paths <- list()

for (d in dirs) {
  sid <- basename(d)
  z <- read_rounds(d)
  if (is.null(z) || nrow(z) == 0) next

  # Provenance guard: rounds spliced from different runs are not comparable.
  span <- as.numeric(difftime(max(z$mtime), min(z$mtime), units = "days"))
  mixed <- span > MAX_MTIME_SPAN_DAYS

  z[, event := classify_event(CHROM, seg_start, seg_end, state)]
  rounds <- sort(unique(z$round))
  final  <- max(rounds)

  ce <- z[!is.na(event)]
  # How many distinct rounds call each event -> corroboration.
  stab <- ce[, .(n_rounds_called = uniqueN(round)), by = event]
  ce   <- merge(ce, stab, by = "event")

  per <- rbindlist(lapply(rounds, function(k) {
    ck <- unique(ce[round == k, .(event, n_rounds_called, LLR)])
    data.table(
      round      = k,
      n_canon    = uniqueN(ck$event),
      n_corrob   = uniqueN(ck[n_rounds_called >= MIN_ROUNDS, event]),
      llr_canon  = sum(ck$LLR, na.rm = TRUE),
      n_segs     = z[round == k & !is.na(state) & state != "neu", .N],
      complete   = all(file.exists(file.path(d, sprintf(REQUIRED_FILES, k)))),
      events     = paste(sort(unique(ck$event)), collapse = ";")
    )
  }))

  # Selection: most corroborated canonical events, then most canonical, then
  # highest canonical LLR, then fewest non-neutral segments (parsimony -- do not
  # reward a round that simply calls more of everything), then earliest round.
  # Completeness first: never hand downstream a round it cannot build.
  setorder(per, -complete, -n_corrob, -n_canon, -llr_canon, n_segs, round)
  sel <- per[1]

  fin <- per[round == final]

  rows[[length(rows) + 1L]] <- data.table(
    sample_id        = sid,
    n_rounds         = length(rounds),
    mixed_provenance = mixed,
    selected_round   = sel$round,
    final_round      = final,
    complete         = sel$complete,
    n_rounds_complete= per[complete == TRUE, .N],
    n_canon_selected = sel$n_canon,
    n_canon_final    = fin$n_canon,
    delta_vs_final   = sel$n_canon - fin$n_canon,
    n_corrob_selected= sel$n_corrob,
    n_segs_selected  = sel$n_segs,
    n_segs_final     = fin$n_segs,
    events_selected  = sel$events,
    events_final     = fin$events,
    n_canon_union    = uniqueN(ce$event),
    per_round_canon  = paste(per[order(round)]$n_canon, collapse = "/"),
    per_round_segs   = paste(per[order(round)]$n_segs,  collapse = "/")
  )

  fp <- sprintf(ROUND_FILES, sel$round)
  paths[[length(paths) + 1L]] <- data.table(
    sample_id = sid, selected_round = sel$round,
    file = fp, path = file.path(d, fp),
    exists = file.exists(file.path(d, fp))
  )
}

out  <- rbindlist(rows)
pout <- rbindlist(paths)
setorder(out, sample_id)

dir.create("results", showWarnings = FALSE)
f1 <- file.path("results", paste0("numbat_selected_round", SUFFIX, ".csv"))
f2 <- file.path("results", paste0("numbat_selected_paths", SUFFIX, ".csv"))
fwrite(out, f1); fwrite(pout, f2)

# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------
clean <- out[mixed_provenance == FALSE]

cat("\n=== per-sample round selection ===\n")
cat("samples:", nrow(out), "  (clean:", nrow(clean),
    " mixed-provenance excluded from totals:", nrow(out) - nrow(clean), ")\n\n")

cat("canonical RB events across the clean cohort:\n")
cat("  numbat final round :", clean[, sum(n_canon_final)], "\n")
cat("  SELECTED round     :", clean[, sum(n_canon_selected)],
    sprintf("  (+%d)", clean[, sum(n_canon_selected) - sum(n_canon_final)]), "\n")
cat("  cross-round union  :", clean[, sum(n_canon_union)],
    "  <- not a numbat state; cannot be built into an object\n")

cat("\nselected round index (watch for collapse onto round 1):\n")
print(clean[, .N, by = selected_round][order(selected_round)])

cat("\nis selection just picking the least-filtered round?\n")
cat("  mean non-neutral segments, selected round:",
    round(clean[, mean(n_segs_selected)], 1), "\n")
cat("  mean non-neutral segments, final round   :",
    round(clean[, mean(n_segs_final)], 1), "\n")

cat("\nsamples where selection beats final:", clean[delta_vs_final > 0, .N],
    " unchanged:", clean[delta_vs_final == 0, .N],
    " worse:", clean[delta_vs_final < 0, .N], "\n")

cat("\nselected round still misses this many union events:",
    clean[, sum(n_canon_union) - sum(n_canon_selected)], "\n")

cat("\nartefact completeness of selected rounds:\n")
cat("  complete   :", out[complete == TRUE, .N], "\n")
cat("  INCOMPLETE :", out[complete == FALSE, .N],
    "(no complete round exists for these)\n")

cat("\nwrote", f1, "\n      ", f2, "\n")
