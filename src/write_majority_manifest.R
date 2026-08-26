# Rewrite results/numbat_selected_round.csv under the majority-of-rounds rule.
#
# Selection rule (replaces the event-maximising rule in select_numbat_round.R):
#   1. An RB event is REPORTED for a sample when a majority of that sample's
#      available consensus rounds call it (arm-coverage-thresholded definition).
#   2. The RDS is then built from the round whose own event set EQUALS that
#      majority set -- so the object always displays exactly what is reported.
#   3. Ties are broken toward numbat's own default i = 2, then the lowest round.
#
# A qualifying round exists for all 36 samples, so no sample needs a spliced
# object (which would be incoherent anyway: tree, clone assignments and per-cell
# posteriors are jointly estimated within one iteration).
#
# The original columns are preserved so existing consumers keep working, but the
# per-round scoring columns (n_canon_selected, events_selected, delta_vs_final,
# ...) were computed under the OLD unthresholded event definition and the old
# rule; they are stale. `selected_round` and the new columns below are current.
suppressPackageStartupMessages(library(data.table))

OUT  <- "results/numbat_selected_round.csv"
BAK  <- "results/numbat_selected_round_premajority_20260826.csv"
stopifnot(file.exists(BAK))

old <- fread(BAK)
new <- fread("results/round_matching_majority.csv")

m <- merge(old, new[, .(sample_id, majority_set, n_majority, best_round,
                        n_rounds_tied, all_best)],
           by = "sample_id", all.x = TRUE)

m[, prev_selected_round := selected_round]
m[!is.na(best_round), selected_round := as.integer(best_round)]
m[, selection_rule := ifelse(is.na(best_round), "unchanged (not scored)",
                             "majority-of-rounds")]
setnames(m, "best_round", "majority_match_round")

# Samples outside the scored set (SRR, and the three held-out SRX) keep whatever
# round they had. Nothing here rebuilds them.
changed <- m[!is.na(majority_match_round) & selected_round != prev_selected_round]

setcolorder(m, c("sample_id", "selected_round", "prev_selected_round",
                 "selection_rule", "majority_set", "n_majority",
                 "majority_match_round", "n_rounds_tied", "all_best"))
fwrite(m, OUT)

cat("wrote", OUT, "  rows:", nrow(m), "\n")
cat("scored (SRX, not held out):", sum(!is.na(m$majority_match_round)), "\n")
cat("round changed vs previous manifest:", nrow(changed), "\n\n")
print(changed[, .(sample_id, prev_selected_round, selected_round, majority_set)])

cat("\nselected_round distribution (scored samples):\n")
print(table(m[!is.na(majority_match_round), selected_round]))

# The rebuild list is defined by disagreement with what is ON DISK, which is the
# previous manifest's round -- not by disagreement with i = 2.
cat("\nsamples needing an RDS rebuild:", nrow(changed), "\n")
writeLines(sort(changed$sample_id), "results/rds_rebuild_list_20260826.txt")
cat("wrote results/rds_rebuild_list_20260826.txt\n")
