# Which samples permit a WITHIN-SAMPLE contrast between clones that differ by an
# RB SCNA of interest (1q+, 2p+, 6p+, 16q-)?
#
# Eligibility here means: two clones exist whose GT_opt token sets are nested
# (A subset of B), and the tokens B adds map to an RB SCNA. That is a
# parent/descendant pair differing by the event of interest, holding the rest of
# the genotype fixed -- the cleanest available internal comparison.
#
# Reads flat per-round TSVs at the SELECTED round; loads no Seurat or numbat
# objects. Cross-checked against config/large_clone_comparisons.yaml.

suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr)})
setwd("/project2/cobrinik_1090/external_rb_scrnaseq_proj")
out_dir <- "results/diploid_audit"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

RB <- c("1q+", "2p+", "6p+", "16q-")
MIN_CELLS <- 20L   # both sides must clear this to be worth a diffex

centromeres_bp <- c(
  "1"=123500000L,"2"=93100000L,"3"=92200000L,"4"=50700000L,"5"=48300000L,
  "6"=59200000L,"7"=59800000L,"8"=45000000L,"9"=44400000L,"10"=40600000L,
  "11"=52700000L,"12"=35600000L,"13"=17000000L,"14"=17100000L,"15"=18400000L,
  "16"=37300000L,"17"=24900000L,"18"=18200000L,"19"=25900000L,"20"=28200000L,
  "21"=11900000L,"22"=14000000L)

scna_label_of <- function(chrom, start, end, cnv_state) {
  chr <- as.character(chrom); mid <- (start + end)/2
  arm <- ifelse(!chr %in% names(centromeres_bp), "?",
                ifelse(mid < centromeres_bp[chr], "p", "q"))
  sfx <- dplyr::case_when(cnv_state %in% c("amp","bamp") ~ "+",
                          cnv_state %in% c("del","bdel") ~ "-",
                          cnv_state %in% c("loh","cnloh") ~ "cnloh",
                          TRUE ~ cnv_state)
  paste0(chr, arm, sfx)
}

sel <- read_csv("results/numbat_selected_round.csv", show_col_types = FALSE) |>
  select(sample_id, selected_round)
samples <- sel$sample_id[grepl("^SRX", sel$sample_id)]
samples <- samples[dir.exists(file.path("output/numbat_sridhar", samples))]

pairs <- list()
for (sid in samples) {
  rnd <- sel$selected_round[sel$sample_id == sid][1]
  cp_f <- sprintf("output/numbat_sridhar/%s/clone_post_%s.tsv", sid, rnd)
  sc_f <- sprintf("output/numbat_sridhar/%s/segs_consensus_%s.tsv", sid, rnd)
  if (!file.exists(cp_f) || !file.exists(sc_f)) next

  cp <- read.delim(cp_f, stringsAsFactors = FALSE, colClasses = "character")
  sc <- read.delim(sc_f, stringsAsFactors = FALSE)
  cp$GT_opt[is.na(cp$GT_opt)] <- ""
  cp$clone_opt <- as.integer(cp$clone_opt)

  sc_nn <- sc[!is.na(sc$cnv_state) & sc$cnv_state != "neu", , drop = FALSE]
  tok_label <- setNames(
    scna_label_of(sc_nn$CHROM, sc_nn$seg_start, sc_nn$seg_end, sc_nn$cnv_state),
    sc_nn$seg_cons)

  cl <- cp |> count(clone_opt, GT_opt, name = "n_cells") |> arrange(clone_opt)
  toks <- lapply(cl$GT_opt, function(g)
    if (g == "") character(0) else str_split_1(g, ","))

  for (i in seq_len(nrow(cl))) for (j in seq_len(nrow(cl))) {
    if (i == j) next
    a <- toks[[i]]; b <- toks[[j]]
    if (!all(a %in% b)) next          # need A nested inside B
    gained <- setdiff(b, a)
    if (!length(gained)) next
    lab <- unique(unname(tok_label[gained]))
    # A clone can gain a NEW SEGMENT on an arm it ALREADY carries -- e.g.
    # SRX10031194 clone 6 -> 7 adds a second 1q token while both clones are
    # already 1q+. Token-level nesting calls that "gains 1q+", but it is not a
    # with/without contrast: both sides carry the event. Require the label to be
    # absent from the preceding clone and present in the descendant.
    lab_a <- unique(unname(tok_label[a]))
    lab_b <- unique(unname(tok_label[b]))
    rb  <- setdiff(intersect(RB, lab_b), lab_a)
    if (!length(rb)) next
    lab <- setdiff(lab_b, lab_a)
    pairs[[length(pairs) + 1L]] <- tibble(
      sample_id = sid, selected_round = rnd,
      clone_wo = cl$clone_opt[i], clone_w = cl$clone_opt[j],
      n_cells_wo = cl$n_cells[i], n_cells_w = cl$n_cells[j],
      wo_is_diploid = length(a) == 0,
      rb_gained = paste(rb, collapse = ";"),
      all_gained = paste(sort(lab[!is.na(lab)]), collapse = ";"),
      background = paste(sort(lab_a[!is.na(lab_a)]), collapse = ";"),
      n_extra_events = length(setdiff(lab, rb)))
  }
}
pr <- bind_rows(pairs)
pr$usable <- pr$n_cells_wo >= MIN_CELLS & pr$n_cells_w >= MIN_CELLS
write_csv(pr, file.path(out_dir, "clone_pairs_rb_scna.csv"))

# Best pair per (sample, RB SCNA): prefer usable, then fewest confounding extra
# events, then most cells on the smaller side.
best <- pr |>
  tidyr::separate_rows(rb_gained, sep = ";") |>
  group_by(sample_id, rb_gained) |>
  arrange(desc(usable), n_extra_events, desc(pmin(n_cells_wo, n_cells_w)),
          .by_group = TRUE) |>
  slice(1) |> ungroup()
write_csv(best, file.path(out_dir, "clone_pairs_rb_scna_best.csv"))

cat("=== samples with a usable within-sample clone pair, by RB SCNA ===\n")
cat("(both clones >=", MIN_CELLS, "cells; n_extra = other events also gained)\n\n")
for (r in RB) {
  b <- best |> filter(rb_gained == r, usable) |> arrange(desc(pmin(n_cells_wo, n_cells_w)))
  cat("--", r, ":", nrow(b), "samples --\n")
  if (nrow(b)) print(as.data.frame(b |> select(sample_id, clone_wo, clone_w,
      n_cells_wo, n_cells_w, wo_is_diploid, n_extra_events)))
  cat("\n")
}
cat("=== union of eligible samples ===\n")
u <- best |> filter(usable) |> group_by(sample_id) |>
  summarise(rb_testable = paste(sort(unique(rb_gained)), collapse = ";"), .groups="drop")
print(as.data.frame(u))
cat("\ntotal eligible samples:", nrow(u), "of", length(samples), "\n")
