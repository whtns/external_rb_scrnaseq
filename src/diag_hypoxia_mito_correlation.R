#!/usr/bin/env Rscript
# Is the hypoxia score inversely related to the mitochondrial score?
#
# Two DIFFERENT questions, only one of which is empirical:
#
#   (a) composite hypoxia_score  vs  mito_score
#       Inverse BY CONSTRUCTION. add_hypoxia_score() computes
#           hypoxia_score = rescale(rowMeans(hypoxia, MT)),  MT stored negated
#       so mito enters the composite with a minus sign. A negative correlation
#       here is arithmetic, not biology.
#
#   (b) raw hypoxia module score  vs  raw mito module score
#       The real question. If these are strongly anti-correlated, the composite
#       is adding two views of ONE signal (and the split's fence is effectively
#       thresholding mitochondrial content twice). If they are independent, the
#       composite genuinely mixes two signals.
#
# Reports both, per sample, at the cell level.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

seus <- sort(Sys.glob("output/seurat/*_seu_hypoxia.rds"))
seus <- seus[!grepl("_labeled", seus)]
cat("samples:", length(seus), "\n\n")

rows <- lapply(seus, function(p) {
  sid <- stringr::str_extract(p, "SR[RX][0-9]+")
  tryCatch({
    m <- readRDS(p)@meta.data
    need <- c("hypoxia", "MT", "hypoxia_score")
    if (!all(need %in% colnames(m))) return(NULL)
    # MT is stored negated; mito_score is its natural orientation.
    mito <- -m$MT
    data.frame(
      sample_id = sid,
      n_cells   = nrow(m),
      # (b) the empirical question
      r_hyp_mito       = suppressWarnings(cor(m$hypoxia, mito, use = "complete.obs")),
      rho_hyp_mito     = suppressWarnings(cor(m$hypoxia, mito, method = "spearman",
                                              use = "complete.obs")),
      # (a) the by-construction one, for contrast
      r_score_mito     = suppressWarnings(cor(m$hypoxia_score, mito, use = "complete.obs")),
      r_score_hyp      = suppressWarnings(cor(m$hypoxia_score, m$hypoxia, use = "complete.obs")),
      # how much of the composite's variance each half carries
      sd_hyp           = sd(m$hypoxia, na.rm = TRUE),
      sd_mito          = sd(mito, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }, error = function(e) { message("!! ", sid, ": ", conditionMessage(e)); NULL })
})

d <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
stopifnot(nrow(d) > 0)

cat("=== (b) EMPIRICAL: raw hypoxia module vs raw mito module, per cell ===\n")
print(d[, c("sample_id", "n_cells", "r_hyp_mito", "rho_hyp_mito")], row.names = FALSE)
cat("\nPearson r  : median ", round(median(d$r_hyp_mito), 3),
    "  range [", round(min(d$r_hyp_mito), 3), ", ",
    round(max(d$r_hyp_mito), 3), "]\n", sep = "")
cat("Spearman rho: median ", round(median(d$rho_hyp_mito), 3),
    "  range [", round(min(d$rho_hyp_mito), 3), ", ",
    round(max(d$rho_hyp_mito), 3), "]\n", sep = "")
cat("samples with r < 0: ", sum(d$r_hyp_mito < 0), "/", nrow(d), "\n", sep = "")

cat("\n=== (a) BY CONSTRUCTION: composite hypoxia_score vs mito / vs hypoxia ===\n")
print(d[, c("sample_id", "r_score_mito", "r_score_hyp")], row.names = FALSE)
cat("\nr(hypoxia_score, mito_score): median ", round(median(d$r_score_mito), 3), "\n", sep = "")
cat("r(hypoxia_score, hypoxia)   : median ", round(median(d$r_score_hyp), 3), "\n", sep = "")

cat("\n=== which half drives the composite? (sd of each ingredient) ===\n")
d$mito_share <- d$sd_mito / (d$sd_hyp + d$sd_mito)
print(d[, c("sample_id", "sd_hyp", "sd_mito", "mito_share")], row.names = FALSE)
cat("\nmito share of ingredient sd: median ",
    round(median(d$mito_share), 3), "  range [",
    round(min(d$mito_share), 3), ", ", round(max(d$mito_share), 3), "]\n", sep = "")

readr::write_csv(d, "results/hypoxia_cluster_split/hypoxia_mito_correlation.csv")
cat("\nwrote results/hypoxia_cluster_split/hypoxia_mito_correlation.csv\n")
cat("DIAG DONE\n")
