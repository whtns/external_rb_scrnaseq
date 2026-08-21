# Confirm each rebuilt *_numbat.rds actually holds its selected consensus round.
# Compares the object's segs_consensus against every segs_consensus_k.tsv on
# disk and reports which round it matches, plus canonical RB event counts.
suppressPackageStartupMessages({library(data.table)})

EXC <- c("SRX11133592", "SRX11133593", "SRX11133594")
man <- fread("results/numbat_selected_round.csv")
man <- man[grepl("^SRX", sample_id) & !sample_id %in% EXC]

canon <- function(sc) {
  if (is.null(sc) || nrow(sc) == 0) return(0L)
  sc <- as.data.table(sc); sc[, CHROM := as.character(CHROM)]
  st <- if ("cnv_state_post" %in% names(sc)) sc$cnv_state_post else sc$cnv_state
  g <- grepl("amp", st); l <- grepl("del|loh", st)
  sum(c(any(sc$CHROM=="1"  & sc$seg_end>125e6  & g, na.rm=TRUE),
        any(sc$CHROM=="2"  & sc$seg_start<93e6 & g, na.rm=TRUE),
        any(sc$CHROM=="6"  & sc$seg_start<59e6 & g, na.rm=TRUE),
        any(sc$CHROM=="13"                     & l, na.rm=TRUE),
        any(sc$CHROM=="16" & sc$seg_end>36.8e6 & l, na.rm=TRUE)))
}

out <- list()
for (j in seq_len(nrow(man))) {
  s <- man$sample_id[j]; k <- as.integer(man$selected_round[j])
  f <- sprintf("output/numbat_sridhar/%s_numbat.rds", s)
  if (!file.exists(f)) next
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) { cat(s, "UNREADABLE\n"); next }
  sc <- tryCatch(as.data.table(nb$segs_consensus), error = function(e) NULL)
  match_k <- NA_integer_
  for (kk in 1:8) {
    p <- sprintf("output/numbat_sridhar/%s/segs_consensus_%d.tsv", s, kk)
    if (!file.exists(p)) next
    d <- fread(p)
    if (!is.null(sc) && nrow(d) == nrow(sc) &&
        isTRUE(all.equal(sort(d$seg_start), sort(sc$seg_start)))) match_k <- kk
  }
  out[[length(out)+1L]] <- data.table(
    sample_id = s, selected = k, rds_matches_round = match_k,
    ok = !is.na(match_k) && match_k == k,
    n_canon = canon(sc),
    mtime = as.character(file.info(f)$mtime))
  rm(nb, sc); invisible(gc(verbose = FALSE))
}
r <- rbindlist(out)
fwrite(r, "results/verify_rebuilt_rds.csv")
cat("\n=== rebuilt RDS verification ===\n")
cat("checked:", nrow(r), "  holding the SELECTED round:", sum(r$ok), "\n")
if (any(!r$ok)) { cat("\nMISMATCHES:\n"); print(r[ok == FALSE]) }
cat("\ncanonical RB events across rebuilt objects:", sum(r$n_canon), "\n")
cat("rebuilt today:", sum(grepl("^2026-08-20", r$mtime)), "of", nrow(r), "\n")
