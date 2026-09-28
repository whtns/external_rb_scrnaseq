#!/usr/bin/env Rscript
# How much does a numbat run's answer depend on WHICH consensus round you read?
#
# The cohort's round is chosen by rule -- round_mode = "final", the highest round
# whose nine artefacts all exist -- and for most samples was never inspected.
# Clone numbering is round-specific, so the choice propagates into every
# downstream clone key. This quantifies, per eligible tumor, how far the rounds
# actually disagree.
#
# Reads only the flat per-round TSVs; no numbat object is built, nothing is
# written outside results/.
#
# Outputs:
#   results/round_variation_by_sample.csv  one row per sample x round
#   results/round_variation_summary.csv    one row per sample
#   results/round_variation_arm_frac.csv   arm coverage per sample x round x event
# Writes: docs/numbat_round_variation.md is the narrative built from these.

setwd("/project2/cobrinik_1090/external_rb_scrnaseq_proj")
suppressMessages(devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers"))

NUMBAT_DIR <- "output/numbat_t1e5"

tri <- read.csv("results/diploid_audit/rb_scna_triage_samples.csv")
el  <- sort(tri$sample_id[tri$best_priority %in% "P1_sufficient"])
act <- read.csv("results/numbat_active_round.csv")
retained <- tri$sample_id[toupper(trimws(as.character(tri$paper_retained))) == "TRUE"]

rows <- list(); af <- list()
for (s in el) {
  d  <- file.path(NUMBAT_DIR, s)
  ks <- numbat_complete_rounds(d)
  a  <- act$round[act$sample_id == s]; a <- if (length(a)) a[1] else NA_integer_
  for (k in ks) {
    segs <- as.data.frame(data.table::fread(
      file.path(d, sprintf("segs_consensus_%d.tsv", k)), showProgress = FALSE))
    cp <- as.data.frame(data.table::fread(
      file.path(d, sprintf("clone_post_%d.tsv", k)), select = "clone_opt", showProgress = FALSE))

    # numbat_rb_arm_frac() needs a data.frame; handing it a path returns all
    # zeros in silence.
    fr <- numbat_rb_arm_frac(segs)
    ev <- sort(names(fr)[fr > 0 & fr >= RB_MIN_ARM_FRAC])

    st <- if ("cnv_state_post" %in% names(segs)) segs$cnv_state_post else segs$cnv_state
    nn <- unique(segs[st != "neu", c("CHROM", "seg_start", "seg_end")])

    rows[[length(rows) + 1L]] <- data.frame(
      sample = s, round = k, active = !is.na(a) && k == a,
      n_clones = length(unique(stats::na.omit(cp$clone_opt))),
      n_events = nrow(nn),
      rb = if (length(ev)) paste(ev, collapse = ";") else "(none)",
      bp = paste(sort(paste0(nn$CHROM, ":", nn$seg_start)), collapse = ","),
      stringsAsFactors = FALSE)
    af[[length(af) + 1L]] <- data.frame(
      sample = s, round = k, active = !is.na(a) && k == a,
      event = names(fr), arm_frac = round(unname(fr), 3),
      called = unname(fr) > 0 & unname(fr) >= RB_MIN_ARM_FRAC,
      stringsAsFactors = FALSE)
  }
}
df <- do.call(rbind, rows)
write.csv(df[, setdiff(names(df), "bp")], "results/round_variation_by_sample.csv", row.names = FALSE)
write.csv(do.call(rbind, af), "results/round_variation_arm_frac.csv", row.names = FALSE)

summ <- do.call(rbind, lapply(split(df, df$sample), function(x) {
  x <- x[order(x$round), ]
  # Breakpoint agreement between CONSECUTIVE rounds, reported as the worst step.
  jac <- NA_real_
  if (nrow(x) > 1) {
    sets <- strsplit(x$bp, ",")
    jac <- min(vapply(2:nrow(x), function(i) {
      a <- sets[[i - 1L]]; b <- sets[[i]]
      if (length(union(a, b)) == 0) 1 else length(intersect(a, b)) / length(union(a, b))
    }, numeric(1)))
  }
  data.frame(
    sample = x$sample[1], rounds = paste(x$round, collapse = ","),
    paper_retained = x$sample[1] %in% retained,
    rb_varies = length(unique(x$rb)) > 1,
    clones_vary = length(unique(x$n_clones)) > 1,
    clone_range = paste0(min(x$n_clones), "-", max(x$n_clones)),
    events_range = paste0(min(x$n_events), "-", max(x$n_events)),
    min_bp_jaccard = round(jac, 3),
    active_round = x$round[x$active][1],
    active_rb = x$rb[x$active][1],
    rb_sets = paste(unique(x$rb), collapse = " | "),
    stringsAsFactors = FALSE)
}))
summ <- summ[order(summ$min_bp_jaccard), ]
write.csv(summ, "results/round_variation_summary.csv", row.names = FALSE)

cat("samples:", nrow(summ), " rb_varies:", sum(summ$rb_varies),
    " clones_vary:", sum(summ$clones_vary),
    " jaccard<0.9:", sum(summ$min_bp_jaccard < 0.9, na.rm = TRUE), "\n\n")
cat("=== samples whose RB call set changes across rounds ===\n")
print(summ[summ$rb_varies, c("sample","rounds","active_round","rb_sets")], row.names = FALSE)
cat("\n=== arm coverage for those samples ===\n")
afd <- do.call(rbind, af)
sel <- afd[afd$sample %in% summ$sample[summ$rb_varies] & afd$arm_frac > 0, ]
print(sel[order(sel$sample, sel$event, sel$round),
          c("sample","round","active","event","arm_frac","called")], row.names = FALSE)
