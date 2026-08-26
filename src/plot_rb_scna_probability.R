# Per-cell SCNA posterior at the canonical RB events, for the samples whose
# event set changed under per-sample round selection (the set in
# results/rb_scna_gain_by_round_selection.csv, i.e. the gained report).
#
# Produces, per sample, the plot_variability_at_SCNA figure -- probability on y,
# cells on x, one facet per supporting segment -- plus a summary row per facet
# so "is this event actually supported per cell" is answerable numerically and
# not only by eye. That distinction matters: segs_consensus can call an event
# whose per-cell posterior sits at zero (SRX10831280's 13q loss does exactly
# this), and only the per-cell view shows it.
#
# Usage: Rscript src/plot_rb_scna_probability.R [manifest_csv] [out_dir]

suppressPackageStartupMessages({
  library(numbat); library(ggplot2); library(patchwork); library(dplyr)
  library(data.table); library(numbatHelpers)
})

args    <- commandArgs(trailingOnly = TRUE)
MAN     <- if (length(args) >= 1) args[[1]] else "results/rb_scna_gain_by_round_selection.csv"
OUT_DIR <- if (length(args) >= 2) args[[2]] else "results/rb_scna_probability"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

samples <- fread(MAN)$sample_id
cat("samples:", length(samples), "\n")

pdfs <- character(0); rows <- list()

for (s in samples) {
  f <- sprintf("output/numbat_sridhar/%s_numbat.rds", s)
  if (!file.exists(f)) { cat("  MISSING", s, "\n"); next }
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) { cat("  UNREADABLE", s, "\n"); next }

  p <- tryCatch(plot_rb_scna_probability(nb, s), error = function(e) {
    cat("  plot failed", s, ":", conditionMessage(e), "\n"); NULL })

  if (!is.null(p)) {
    out <- file.path(OUT_DIR, paste0(s, "_rb_scna_probability.pdf"))
    ok <- tryCatch({ ggsave(out, p, width = 14, height = 9, limitsize = FALSE); TRUE },
                   error = function(e) { cat("  save failed", s, ":", conditionMessage(e), "\n"); FALSE })
    if (ok) { pdfs <- c(pdfs, out); cat(sprintf("  %-14s %d facets\n", s, nlevels(p$data$seg))) }

    # Quantify each facet: how much of the tumour actually carries the call.
    d <- as.data.table(p$data)
    rows[[length(rows) + 1L]] <- d[, .(
      sample_id   = s,
      n_cells     = uniqueN(cell),
      frac_p_gt90 = round(mean(p_cnv > 0.9, na.rm = TRUE), 3),
      frac_p_gt50 = round(mean(p_cnv > 0.5, na.rm = TRUE), 3),
      median_p    = round(stats::median(p_cnv, na.rm = TRUE), 3)
    ), by = .(facet = as.character(seg), event, cnv_state)]
  }
  rm(nb); invisible(gc(verbose = FALSE))
}

if (length(rows)) {
  tab <- rbindlist(rows)
  setorder(tab, sample_id, event)
  fwrite(tab, file.path("results", "rb_scna_probability_summary.csv"))
  cat("\nwrote results/rb_scna_probability_summary.csv -", nrow(tab), "facets\n")
  cat("\nfacets whose call is NOT supported per cell (frac_p_gt90 < 0.05):\n")
  print(tab[frac_p_gt90 < 0.05], row.names = FALSE)
}

if (length(pdfs)) {
  comb <- file.path("results", "rb_scna_probability_gained.pdf")
  qpdf::pdf_combine(pdfs, output = comb)
  cat("\nwrote", comb, "-", qpdf::pdf_length(comb), "pages\n")
}
