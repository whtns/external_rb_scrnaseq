# Verify each *_numbat.rds on disk holds the round the manifest now specifies.
#
# Independent of rebuild history: reads the object, compares its segs_consensus
# breakpoints against every segs_consensus_k.tsv on disk, and reports which
# round it actually matches. Also checks the object carries a usable phylogeny,
# since building from a round whose tree artefacts are absent yields an object
# that loads fine but has no tree.
#
# Writes results/verify_majority_rds.csv and the list of samples still needing
# a rebuild to results/rds_rebuild_list_pass2.txt.
suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

NB_DIR <- "output/numbat_sridhar"
EXC    <- c("SRX11133592", "SRX11133593", "SRX11133594")

man <- fread("results/numbat_selected_round.csv")
man <- man[grepl("^SRX", sample_id) & !sample_id %in% EXC]

out <- list()
for (j in seq_len(nrow(man))) {
  s <- man$sample_id[j]
  k <- as.integer(man$selected_round[j])
  f <- file.path(NB_DIR, sprintf("%s_numbat.rds", s))
  if (!file.exists(f)) { cat(s, "MISSING RDS\n"); next }
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) { cat(s, "UNREADABLE\n"); next }

  sc <- tryCatch(as.data.table(nb$segs_consensus), error = function(e) NULL)
  match_k <- NA_integer_
  for (kk in 1:12) {
    p <- file.path(NB_DIR, s, sprintf("segs_consensus_%d.tsv", kk))
    if (!file.exists(p)) next
    d <- fread(p)
    if (!is.null(sc) && nrow(d) == nrow(sc) &&
        isTRUE(all.equal(sort(d$seg_start), sort(sc$seg_start)))) match_k <- kk
  }

  has_tree  <- !is.null(nb$gtree) || !is.null(nb$treeML)
  has_clone <- !is.null(nb$clone_post) && NROW(nb$clone_post) > 0
  ev <- tryCatch(numbat_rb_events(nb$segs_consensus), error = function(e) character(0))

  out[[length(out) + 1L]] <- data.table(
    sample_id     = s,
    manifest_round = k,
    rds_round      = match_k,
    round_ok       = !is.na(match_k) && match_k == k,
    has_tree       = has_tree,
    has_clone_post = has_clone,
    n_rb           = length(ev),
    rb_events      = paste(sort(ev), collapse = ";"),
    majority_set   = man$majority_set[j],
    events_match   = identical(sort(ev),
                     sort(if (nzchar(man$majority_set[j]))
                            strsplit(man$majority_set[j], ";")[[1]] else character(0))))
  rm(nb, sc); invisible(gc(verbose = FALSE))
}
R <- rbindlist(out)
fwrite(R, "results/verify_majority_rds.csv")

cat("\n=== verification ===\n")
cat("checked:", nrow(R), "\n")
cat("  RDS holds the manifest round      :", sum(R$round_ok), "\n")
cat("  object carries a phylogeny        :", sum(R$has_tree), "\n")
cat("  object carries clone_post         :", sum(R$has_clone_post), "\n")
cat("  event set matches the majority set:", sum(R$events_match), "\n")

bad <- R[round_ok == FALSE | has_tree == FALSE | has_clone_post == FALSE |
         events_match == FALSE]
if (nrow(bad)) {
  cat("\nNEEDS ATTENTION:\n")
  print(bad[, .(sample_id, manifest_round, rds_round, has_tree, has_clone_post,
                events_match, rb_events, majority_set)])
  writeLines(sort(unique(bad$sample_id)), "results/rds_rebuild_list_pass2.txt")
  cat("\nwrote results/rds_rebuild_list_pass2.txt (", nrow(bad), "samples )\n")
} else {
  cat("\nall objects consistent with the manifest.\n")
  if (file.exists("results/rds_rebuild_list_pass2.txt"))
    file.remove("results/rds_rebuild_list_pass2.txt")
}
cat("\ntotal canonical RB events across the cohort:", sum(R$n_rb), "\n")
