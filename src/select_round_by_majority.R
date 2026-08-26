# Choose the consensus round to build each sample's *_numbat.rds from.
#
# Replaces the event-maximising rule in select_numbat_round.R. See
# docs/numbat_round_policy.md sections 4b/4c for the evidence.
#
# Rule, in precedence order:
#   1. An RB event is REPORTED when a majority of the sample's available
#      consensus rounds call it (arm-coverage-thresholded definition).
#   2. Prefer the round whose own event set best matches that majority set.
#   3. Among rounds tied on (2), prefer one whose full artefact set is on disk.
#      numbat does not always write the tree artefacts for every round --
#      SRX10031191 and SRX11133586 have a complete round 1 and a round 2 that
#      is missing clone_post/treeML/mut_graph/tree_final -- and building from an
#      incomplete round yields an object with no phylogeny, which silently
#      breaks clone trees, heatmaps and clone-based diffex.
#   4. Among rounds still tied, prefer numbat's own default i = 2, then lowest.
#
# Writes results/round_matching_majority.csv.
suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

NB_DIR <- "output/numbat_sridhar"
EXC    <- c("SRX11133592", "SRX11133593", "SRX11133594")

# the nine files process_numbat_rds.R requires for a given round
round_files <- function(k) {
  c(sprintf("segs_consensus_%d.tsv", k), sprintf("clone_post_%d.tsv", k),
    sprintf("joint_post_%d.tsv", k), sprintf("exp_post_%d.tsv", k),
    sprintf("allele_post_%d.tsv", k), sprintf("geno_%d.tsv", k),
    sprintf("treeML_%d.rds", k), sprintf("mut_graph_%d.rds", k),
    sprintf("tree_final_%d.rds", k))
}
is_complete <- function(s, k) all(file.exists(file.path(NB_DIR, s, round_files(k))))

BY <- fread("results/round_convergence_by_round.csv")
ev <- function(s, k) {
  r <- BY[sample_id == s & round == k]
  if (!nrow(r) || !nzchar(r$rb_events)) return(character(0))
  strsplit(r$rb_events, ";")[[1]]
}

samples <- sort(unique(BY$sample_id))
samples <- samples[grepl("^SRX", samples) & !samples %in% EXC]

out <- list()
for (s in samples) {
  ks <- sort(BY[sample_id == s, round])
  n  <- length(ks)
  tally <- table(unlist(lapply(ks, function(k) ev(s, k))))
  maj   <- sort(names(tally)[tally > n / 2])

  disagree <- vapply(ks, function(k) {
    e <- sort(ev(s, k)); length(setdiff(maj, e)) + length(setdiff(e, maj))
  }, integer(1))
  best <- ks[disagree == min(disagree)]

  complete <- vapply(best, function(k) is_complete(s, k), logical(1))
  # step 3: restrict to complete rounds if any tied round is complete
  cand <- if (any(complete)) best[complete] else best
  # step 4
  pick <- if (2L %in% cand) 2L else min(cand)

  out[[length(out) + 1L]] <- data.table(
    sample_id      = s,
    majority_set   = paste(maj, collapse = ";"),
    n_majority     = length(maj),
    best_round     = pick,
    disagreement   = min(disagree),
    exact_match    = min(disagree) == 0L,
    round_complete = is_complete(s, pick),
    n_rounds_tied  = length(best),
    all_best       = paste(best, collapse = ","),
    tied_complete  = paste(best[complete], collapse = ","),
    dropped_for_artefacts = any(!complete) && !identical(cand, best))
}
R <- rbindlist(out)
fwrite(R, "results/round_matching_majority.csv")

cat("samples:", nrow(R), "\n")
cat("  a round matches the majority set exactly:", sum(R$exact_match), "\n")
cat("  chosen round has a complete artefact set:", sum(R$round_complete), "\n")
cat("  choice moved off a tied round because its artefacts were incomplete:",
    sum(R$dropped_for_artefacts), "\n\n")
if (any(R$dropped_for_artefacts)) {
  cat("artefact-driven choices:\n")
  print(R[dropped_for_artefacts == TRUE,
          .(sample_id, majority_set, all_best, tied_complete, best_round)])
}
if (any(!R$round_complete)) {
  cat("\nWARNING - chosen round is INCOMPLETE for:\n")
  print(R[round_complete == FALSE, .(sample_id, best_round, all_best)])
}
cat("\nchosen round distribution:\n"); print(table(R$best_round))
cat("\nwrote results/round_matching_majority.csv\n")
