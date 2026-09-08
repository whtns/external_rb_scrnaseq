# Audit: which tumors have an initial clone that already carries an SCNA,
# and what that means for the consensus diploid object (github #44).
#
# Why this matters. assemble_diploid_seu() pools cells where
#   is.na(seu$scna) | seu$scna == ""
# and seu$scna is derived, per cell, from that cell's numbat GT_opt
# (plot_functions_3.R:150-171). So a cell reaches the diploid pool by exactly
# three routes:
#
#   A. its clone has an empty GT_opt          -> numbat's normal clone. Intended.
#   B. it is absent from clone_post           -> scna is NA. Silent.
#   C. its clone has a NON-empty GT_opt whose tokens do not resolve against the
#      clone-simplification key -> simplify_gt_col() falls back to "" when the
#      join drops every token (plot_functions_15.R:246-247), which is the SAME
#      value a genuinely diploid clone gets. Silent, and aneuploid.
#
# Route A is the only one that means "diploid". B and C put tumor cells into a
# reference meant to be copy-neutral, so this script counts all three per sample.
#
# Reads flat per-round TSVs, never the GB-scale *_numbat.rds, so it is cheap.

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(DBI); library(readr)
})

proj <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj"
setwd(proj)
out_dir <- "results/diploid_audit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# hg38 centromere midpoints, copied from compute_clone_simplifications()
# (metadata_functions_1.R:474) so the arm labels here match the pipeline's.
centromeres_bp <- c(
  "1" = 123500000L, "2" = 93100000L,  "3" = 92200000L,  "4" = 50700000L,
  "5" = 48300000L,  "6" = 59200000L,  "7" = 59800000L,  "8" = 45000000L,
  "9" = 44400000L,  "10"= 40600000L,  "11"= 52700000L,  "12"= 35600000L,
  "13"= 17000000L,  "14"= 17100000L,  "15"= 18400000L,  "16"= 37300000L,
  "17"= 24900000L,  "18"= 18200000L,  "19"= 25900000L,  "20"= 28200000L,
  "21"= 11900000L,  "22"= 14000000L
)

scna_label_of <- function(chrom, start, end, cnv_state) {
  chr <- as.character(chrom)
  mid <- (start + end) / 2
  arm <- ifelse(!chr %in% names(centromeres_bp), "?",
                ifelse(mid < centromeres_bp[chr], "p", "q"))
  suffix <- dplyr::case_when(
    cnv_state %in% c("amp", "bamp")  ~ "+",
    cnv_state %in% c("del", "bdel")  ~ "-",
    cnv_state %in% c("loh", "cnloh") ~ "cnloh",
    TRUE ~ cnv_state
  )
  paste0(chr, arm, suffix)
}

# ---- inputs -----------------------------------------------------------------

sel <- read_csv("results/numbat_selected_round.csv", show_col_types = FALSE) |>
  select(sample_id, selected_round)

manual_yaml <- yaml::read_yaml("config/large_clone_simplifications.yaml")

panel <- read_csv("data/diploid_panel_samples.csv", show_col_types = FALSE) |>
  select(sample_id, include, panel_reason = reason)

con <- dbConnect(RSQLite::SQLite(), "batch_hashes.sqlite")
low <- dbGetQuery(con, "select filepath, sample_id, cells from seu_cells
                        where seu_type = 'hypoxia_low'")
dbDisconnect(con)
# the _recompute_test copy is a scratch duplicate, not a pipeline input
low <- low[!grepl("_recompute_test", low$filepath), ]
low_cells <- setNames(strsplit(low$cells, "\n", fixed = TRUE), low$sample_id)

samples <- sel$sample_id[grepl("^SRX", sel$sample_id)]
samples <- samples[dir.exists(file.path("output/numbat_sridhar", samples))]
cat("SRX samples with a numbat directory:", length(samples), "\n\n")

# ---- per-sample audit -------------------------------------------------------

clone_rows  <- list()
sample_rows <- list()

for (sid in samples) {
  rnd <- sel$selected_round[sel$sample_id == sid][1]
  cp_f <- sprintf("output/numbat_sridhar/%s/clone_post_%s.tsv", sid, rnd)
  sc_f <- sprintf("output/numbat_sridhar/%s/segs_consensus_%s.tsv", sid, rnd)
  if (!file.exists(cp_f) || !file.exists(sc_f)) {
    cat(sid, ": missing round", rnd, "TSVs -- skipped\n"); next
  }

  cp <- read.delim(cp_f, stringsAsFactors = FALSE, colClasses = "character")
  sc <- read.delim(sc_f, stringsAsFactors = FALSE)
  cp$GT_opt[is.na(cp$GT_opt)] <- ""
  cp$clone_opt <- as.integer(cp$clone_opt)

  # token -> readable event, keyed on seg_cons (GT_opt tokens ARE seg_cons,
  # not seg; seg is shared by the amp/del alternatives of one physical segment)
  sc_nn <- sc[!is.na(sc$cnv_state) & sc$cnv_state != "neu", , drop = FALSE]
  tok_label <- setNames(
    scna_label_of(sc_nn$CHROM, sc_nn$seg_start, sc_nn$seg_end, sc_nn$cnv_state),
    sc_nn$seg_cons
  )

  # Reproduce the key simplify_gt_col() actually joins against: one
  # representative SEG per scna_label (compute_clone_simplifications), joined
  # by "seg". Tokens are seg_cons, so a token resolves only when it also
  # appears as a seg value in the key.
  # Manual YAML entries override computed ones by LABEL (modifyList), so a
  # curated "16q-: 16b" REPLACES the computed representative for 16q- rather
  # than adding to it -- the key can get narrower, not just wider. Reproduce
  # the merge exactly instead of assuming.
  computed_key <- list()
  if (nrow(sc_nn) > 0) {
    k <- sc_nn
    k$scna_label <- scna_label_of(k$CHROM, k$seg_start, k$seg_end, k$cnv_state)
    k <- k[order(k$scna_label, k$seg), ]
    k <- k[!duplicated(k$scna_label), ]
    computed_key <- as.list(as.character(k$seg))
    names(computed_key) <- k$scna_label
  }
  manual_key <- manual_yaml[[sid]]
  merged_key <- if (is.null(manual_key) || length(manual_key) == 0) computed_key
                else modifyList(computed_key, manual_key)
  key_seg <- unlist(merged_key, use.names = FALSE)
  key_seg <- as.character(key_seg[!is.na(key_seg) & key_seg != ""])

  clones <- cp |>
    count(clone_opt, GT_opt, name = "n_cells") |>
    arrange(clone_opt)

  compartment <- cp |>
    group_by(clone_opt) |>
    summarise(compartment = names(sort(table(compartment_opt), decreasing = TRUE))[1],
              .groups = "drop")
  clones <- left_join(clones, compartment, by = "clone_opt")

  clones$n_tokens <- ifelse(clones$GT_opt == "", 0L,
                            lengths(strsplit(clones$GT_opt, ",", fixed = TRUE)))
  clones$events <- vapply(clones$GT_opt, function(g) {
    if (g == "") return("")
    toks <- strsplit(g, ",", fixed = TRUE)[[1]]
    lab <- tok_label[toks]
    paste(sort(unique(ifelse(is.na(lab), paste0("?", toks), lab))), collapse = ";")
  }, character(1))

  # route C: non-empty GT that simplify_gt_col() would collapse to ""
  clones$n_tokens_resolved <- vapply(clones$GT_opt, function(g) {
    if (g == "") return(0L)
    sum(strsplit(g, ",", fixed = TRUE)[[1]] %in% key_seg)
  }, integer(1))
  clones$mislabeled_diploid <- clones$GT_opt != "" & clones$n_tokens_resolved == 0L

  clones$sample_id <- sid
  clones$selected_round <- rnd
  clone_rows[[sid]] <- clones

  # ---- diploid-pool composition for this sample -----------------------------
  norm_bc <- function(x) sub(".", "-", x, fixed = TRUE)
  cp_cell <- norm_bc(cp$cell)
  gt_by_cell <- setNames(cp$GT_opt, cp_cell)
  mis_gt <- clones$GT_opt[clones$mislabeled_diploid]

  lc <- low_cells[[sid]]
  if (is.null(lc)) {
    n_low <- NA_integer_; nA <- NA_integer_; nB <- NA_integer_; nC <- NA_integer_
  } else {
    lc <- norm_bc(lc)
    n_low <- length(lc)
    gt <- gt_by_cell[lc]
    nB <- sum(is.na(gt))                       # absent from clone_post
    nA <- sum(!is.na(gt) & gt == "")           # true diploid clone
    nC <- sum(!is.na(gt) & gt != "" & gt %in% mis_gt)
  }

  clone1_gt <- clones$GT_opt[clones$clone_opt == 1]
  clone1_gt <- if (length(clone1_gt) == 0) NA_character_ else clone1_gt[1]
  tumor <- clones[clones$GT_opt != "", , drop = FALSE]
  founder <- if (nrow(tumor) == 0) NA_character_ else
    tumor$events[which.min(tumor$n_tokens)][1]

  sample_rows[[sid]] <- tibble(
    sample_id = sid, selected_round = rnd,
    n_clones = nrow(clones), n_tumor_clones = nrow(tumor),
    n_cells_clone_post = nrow(cp),
    clone1_gt = clone1_gt,
    initial_clone_has_scna = !is.na(clone1_gt) & clone1_gt != "",
    has_diploid_clone = any(clones$GT_opt == ""),
    n_cells_diploid_clone = sum(clones$n_cells[clones$GT_opt == ""]),
    founder_events = founder,
    n_clones_mislabeled_diploid = sum(clones$mislabeled_diploid),
    n_cells_mislabeled_diploid = sum(clones$n_cells[clones$mislabeled_diploid]),
    lowhyp_cells = n_low,
    pool_A_true_diploid = nA, pool_B_absent_na = nB, pool_C_mislabeled = nC
  )

  cat(sprintf("%-12s round %s  clones %2d (tumor %2d)  clone1 %-8s  diploid clone %-5s  mislabeled clones %d\n",
              sid, rnd, nrow(clones), nrow(tumor),
              ifelse(is.na(clone1_gt), "NA", ifelse(clone1_gt == "", "(empty)", "SCNA")),
              any(clones$GT_opt == ""), sum(clones$mislabeled_diploid)))
}

clone_tbl  <- bind_rows(clone_rows) |>
  select(sample_id, selected_round, clone_opt, GT_opt, events, n_cells,
         compartment, n_tokens, n_tokens_resolved, mislabeled_diploid)
sample_tbl <- bind_rows(sample_rows) |> left_join(panel, by = "sample_id")

sample_tbl$pool_total <- sample_tbl$pool_A_true_diploid +
  sample_tbl$pool_B_absent_na + sample_tbl$pool_C_mislabeled
sample_tbl$pool_pct_contaminated <- round(
  100 * (sample_tbl$pool_B_absent_na + sample_tbl$pool_C_mislabeled) /
    pmax(sample_tbl$pool_total, 1), 1)

write_csv(clone_tbl,  file.path(out_dir, "initial_clone_scna_clones.csv"))
write_csv(sample_tbl, file.path(out_dir, "initial_clone_scna_audit.csv"))

# ---- summary ----------------------------------------------------------------

cat("\n================ RESULT ================\n")
cat("samples audited:", nrow(sample_tbl), "\n\n")

cat("-- 1. does the initial clone (clone_opt == 1) carry an SCNA? --\n")
print(sample_tbl |> count(initial_clone_has_scna))
bad1 <- sample_tbl |> filter(initial_clone_has_scna)
if (nrow(bad1)) print(bad1 |> select(sample_id, include, clone1_gt, founder_events,
                                     n_cells_diploid_clone) |> as.data.frame())

cat("\n-- 2. is there ANY diploid (empty-GT) clone at all? --\n")
print(sample_tbl |> count(has_diploid_clone))
bad2 <- sample_tbl |> filter(!has_diploid_clone)
if (nrow(bad2)) print(bad2 |> select(sample_id, include, n_tumor_clones,
                                     founder_events, lowhyp_cells) |> as.data.frame())

cat("\n-- 3. clones a non-empty GT would still label \"\" (route C) --\n")
bad3 <- sample_tbl |> filter(n_clones_mislabeled_diploid > 0)
cat("samples affected:", nrow(bad3), "of", nrow(sample_tbl), "\n")
if (nrow(bad3)) print(bad3 |> select(sample_id, include, n_clones_mislabeled_diploid,
                                     n_cells_mislabeled_diploid) |> as.data.frame())

cat("\n-- 4. diploid pool composition, panel samples only (include == TRUE) --\n")
inc <- sample_tbl |> filter(!is.na(include), as.logical(include))
print(inc |> select(sample_id, lowhyp_cells, pool_A_true_diploid,
                    pool_B_absent_na, pool_C_mislabeled, pool_total,
                    pool_pct_contaminated) |> as.data.frame())
cat("\ncohort pool: A(true diploid)", sum(inc$pool_A_true_diploid, na.rm = TRUE),
    " B(absent/NA)", sum(inc$pool_B_absent_na, na.rm = TRUE),
    " C(mislabeled)", sum(inc$pool_C_mislabeled, na.rm = TRUE), "\n")

cat("\n-- 5. samples in seus_low_hypoxia but absent from the panel CSV --\n")
missing <- setdiff(names(low_cells), panel$sample_id)
cat(if (length(missing)) paste(missing, collapse = ", ") else "(none)", "\n")
cat("\nwrote", file.path(out_dir, "initial_clone_scna_audit.csv"), "\n")
