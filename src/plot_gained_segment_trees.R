# Segment-labeled clone trees for the samples whose canonical RB event set
# changed under per-sample round selection (results/rb_scna_gain_by_round_selection.csv).
#
# Uses the established path: plot_clone_tree(..., clone_simplifications = NULL),
# which is what the *_segment_tree.pdf targets emit -- edges carry numbat's RAW
# segment labels rather than simplified SCNA names.
#
# Driven from clone_post rather than a Seurat object. save_clone_tree_from_path()
# reads a Seurat file only to restrict to retained cells; with no filtering the
# tree is the unfiltered one, and that keeps this runnable for any sample with a
# numbat RDS regardless of which Seurat objects exist.
#
# A caption maps the cryptic segment labels (1e, 16g, 13b) onto the canonical RB
# arms they support, so the tree can be read without cross-referencing
# segs_consensus by hand.
#
# Usage: Rscript src/plot_gained_segment_trees.R [manifest_csv] [out_dir]

suppressPackageStartupMessages({
  library(numbat); library(ggplot2); library(dplyr); library(tidygraph)
  library(igraph); library(ggtree); library(ggraph)
  library(data.table); library(numbatHelpers)
})

args    <- commandArgs(trailingOnly = TRUE)
MAN     <- if (length(args) >= 1) args[[1]] else "results/rb_scna_gain_by_round_selection.csv"
OUT_DIR <- if (length(args) >= 2) args[[2]] else "results/rb_segment_trees"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

samples <- fread(MAN)$sample_id
cat("samples:", length(samples), "\n")

pdfs <- character(0); notes <- list()

for (s in samples) {
  f <- sprintf("output/numbat_sridhar/%s_numbat.rds", s)
  if (!file.exists(f)) { cat("  MISSING", s, "\n"); next }

  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) { cat("  UNREADABLE", s, "\n"); next }

  cd <- tryCatch(
    nb$clone_post %>% dplyr::select(cell, clone_opt) %>% dplyr::distinct(),
    error = function(e) NULL)
  if (is.null(cd) || nrow(cd) == 0) {
    cat("  no clone_post", s, "\n"); rm(nb); invisible(gc(verbose = FALSE)); next
  }

  # Which raw segments carry a canonical RB event, so the caption can decode them.
  rb  <- tryCatch(numbat_rb_segments(nb$segs_consensus), error = function(e) NULL)
  key <- if (!is.null(rb) && nrow(rb)) {
    paste(vapply(seq_len(nrow(rb)), function(i)
      sprintf("%s = %s%s", rb$seg[i], sub("_.*", "", rb$event[i]),
              ifelse(grepl("gain$", rb$event[i]), "+", "-")), character(1)),
      collapse = "   ")
  } else "no canonical RB event above the arm-coverage floor"

  # Wrap: samples with many supporting segments (SRX11133587 has 10) produce a
  # caption wider than the page, and ggplot silently clips it rather than
  # wrapping, so the last entries vanish off the right edge.
  key_wrapped <- paste(strwrap(key, width = 62), collapse = "\n")

  p <- tryCatch(
    plot_clone_tree(cd, tumor_id = s, nb_path = f, clone_simplifications = NULL,
                    sample_id = s, legend = FALSE, horizontal = FALSE),
    error = function(e) { cat("  tree failed", s, ":", conditionMessage(e), "\n"); NULL })

  if (!is.null(p)) {
    p <- p + ggplot2::labs(
      caption = paste0("canonical RB segments:\n", key_wrapped,
                       "\ncells: ", nrow(cd), "   clones: ", length(unique(cd$clone_opt)))) +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, hjust = 0))
    out <- file.path(OUT_DIR, paste0(s, "_segment_tree.pdf"))
    ok <- tryCatch({ ggsave(out, p, width = 7, height = 7.5); TRUE },
                   error = function(e) { cat("  save failed", s, ":", conditionMessage(e), "\n"); FALSE })
    if (ok) {
      pdfs <- c(pdfs, out)
      cat(sprintf("  %-14s %d clones   %s\n", s, length(unique(cd$clone_opt)), key))
      notes[[length(notes) + 1L]] <- data.table(
        sample_id = s, n_cells = nrow(cd),
        n_clones = length(unique(cd$clone_opt)), rb_segment_key = key)
    }
  }
  rm(nb); invisible(gc(verbose = FALSE))
}

if (length(notes)) {
  fwrite(rbindlist(notes), "results/rb_segment_tree_key.csv")
  cat("\nwrote results/rb_segment_tree_key.csv\n")
}
if (length(pdfs)) {
  comb <- "results/rb_segment_trees_gained.pdf"
  qpdf::pdf_combine(pdfs, output = comb)
  cat("wrote", comb, "-", qpdf::pdf_length(comb), "pages\n")
}
