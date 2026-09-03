#!/usr/bin/env Rscript
# Smoke test for plot_numbat_bulk_clones() (github #41 / #42).
#
# Picks the two samples with the largest known round mismatch plus one control,
# renders the panel, and checks the thing that actually matters: that the round
# the object holds -- and therefore the events the panel draws -- is
# `selected_round` from the manifest, NOT `final_round` (which is what the old
# bulk_clones_final.png path was showing).
#
# Writes only into results/numbat_bulk_clones_SMOKE/. Touches no targets store.

suppressPackageStartupMessages({
  library(devtools)
})
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)

SAMPLES  <- c("SRX10264524", "SRX22868104", "SRX10264519")
OUT_DIR  <- "results/numbat_bulk_clones_SMOKE"
MANIFEST <- "results/numbat_selected_round.csv"

man <- as.data.frame(data.table::fread(MANIFEST))

cat("=== rendering ===\n")
rows <- lapply(SAMPLES, function(s) {
  rds <- file.path("output/numbat_sridhar", paste0(s, "_numbat.rds"))
  if (!file.exists(rds)) { cat("  ", s, ": no RDS\n"); return(NULL) }

  t0  <- Sys.time()
  pdf <- plot_numbat_bulk_clones(rds, out_dir = OUT_DIR, manifest = MANIFEST)
  el  <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  nb  <- readRDS(rds)
  fld <- function(n) tryCatch(nb[[n]], error = function(e) NULL)

  held <- fld("selected_round")
  held <- if (is.null(held)) numbat_match_round(fld("segs_consensus"),
                                                sub("_numbat\\.rds$", "", rds))
          else as.integer(held)

  hit <- man[man$sample_id == s, , drop = FALSE]
  ev_obj <- paste(sort(numbat_rb_events(fld("segs_consensus"))), collapse = ";")

  cat(sprintf("  %s: %s (%ss)\n", s,
              if (is.na(pdf)) "FAILED" else basename(pdf), el))

  data.frame(
    sample_id      = s,
    pdf_ok         = !is.na(pdf) && file.exists(pdf),
    n_pages        = if (!is.na(pdf) && file.exists(pdf))
                       tryCatch(qpdf::pdf_length(pdf), error = function(e) NA_integer_)
                     else NA_integer_,
    round_held     = held,
    round_selected = if (nrow(hit)) as.integer(hit$selected_round[1]) else NA_integer_,
    round_final    = if (nrow(hit)) as.integer(hit$final_round[1])    else NA_integer_,
    events_object  = ev_obj,
    # majority_set, NOT events_selected. events_selected disagrees with the
    # selected round's own events for 25/36 samples -- it is a stale column.
    # majority_set matches the per-round table for 36/36.
    events_sel     = if (nrow(hit)) paste(sort(strsplit(hit$majority_set[1], ";")[[1]]), collapse = ";") else NA,
    events_final   = if (nrow(hit)) paste(sort(strsplit(hit$events_final[1], ";")[[1]]), collapse = ";") else NA,
    stringsAsFactors = FALSE
  )
})

d <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

cat("\n=== reconciliation ===\n")
print(d[, c("sample_id", "pdf_ok", "n_pages",
            "round_held", "round_selected", "round_final")], row.names = FALSE)

cat("\n=== events: object vs manifest ===\n")
for (i in seq_len(nrow(d))) {
  cat(sprintf("%s\n  object   : %s\n  selected : %s\n  final    : %s\n",
              d$sample_id[i], d$events_object[i], d$events_sel[i], d$events_final[i]))
}

ok_round  <- all(d$round_held == d$round_selected, na.rm = TRUE)
ok_pdf    <- all(d$pdf_ok)
ok_pages  <- isTRUE(all(!is.na(d$n_pages)) && all(d$n_pages == 2))
not_final <- any(d$round_selected != d$round_final, na.rm = TRUE)
ev_match  <- d$events_object == d$events_sel

cat("\n=== verdict ===\n")
cat("all PDFs rendered            : ", ok_pdf, "\n", sep = "")
cat("two pages each               : ", ok_pages, "\n", sep = "")
cat("round held == selected_round : ", ok_round, "\n", sep = "")
cat("selected != final round      : ", not_final,
    "  (so this is a real change)\n", sep = "")

ok_events <- all(ev_match, na.rm = TRUE)
cat("events match majority_set    : ", ok_events, "\n", sep = "")
if (!ok_events) {
  cat("\nobject events differ from majority_set for: ",
      paste(d$sample_id[!ev_match], collapse = ", "), "\n", sep = "")
}

if (!(ok_pdf && ok_pages && ok_round && ok_events)) {
  cat("\nSMOKE FAILED\n"); quit(status = 1)
}
cat("\nSMOKE PASSED\n")
