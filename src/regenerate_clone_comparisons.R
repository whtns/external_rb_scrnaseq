# Regenerate config/large_clone_comparisons.yaml from the current numbat objects.
#
# WHY. The curated file was written against pre-revamp numbat objects. The
# majority-of-rounds selection changed the selected round for many samples, and
# both clone numbering and segment labels are round-specific. Validating all 53
# curated entries against the current objects: 0 match exactly, 36 mismatch, 10
# partial, 4 name a clone that no longer exists. e.g. SRX10264523 "5_v_4_2p+" --
# in round 2 clones 5 vs 4 differ by 19d (chr19); the 2p acquisition is 4 -> 7.
#
# CONTRACT the consumers rely on (numbatHelpers::plot_scna_two_clone_res_collages,
# and the str_detect(names(...), scna_of_interest) idiom used throughout):
#   key    "<B>_v_<A>_<tok>[_<tok>...]"  B = acquiring clone, A = its parent
#   token  signed SCNA, e.g. "1q+", "2p+", "6p+", "16q-"
#   value  the segment label(s) supporting those tokens
# The SCNA is matched as a fixed substring of the key, so a transition acquiring
# two events is named once and appears in both SCNAs' collages -- the existing
# "2_v_1_1q+_16q-" convention, preserved here.
#
# NOT REGENERATED, preserved verbatim from the existing file:
#   - SRX11133592/93/94, held out per AGENTS.md and absent from the majority
#     table; their numbat output is untouched so their entries still stand.
#   - "<sample>_branch_<N>" keys. Those address cell-subset objects which need
#     not contain every clone, so a base-sample edge cannot be assumed valid for
#     them. They remain STALE and are reported as such -- review separately.
#
# Writes config/large_clone_comparisons_regenerated.yaml. Does not overwrite the
# live config; diff first.
#
# Usage: Rscript src/regenerate_clone_comparisons.R

suppressPackageStartupMessages(library(data.table))

NB        <- "output/numbat_sridhar"
MIN_CELLS <- 20
MIN_ARM   <- 0.15   # RB_MIN_ARM_FRAC
OUT       <- "config/large_clone_comparisons_regenerated.yaml"

# hg38 arm windows, matching numbatHelpers::.RB_WINDOWS
EV <- data.table(
  token = c("1q+",  "2p+", "6p+", "13q-", "16q-"),
  CHROM = c("1",    "2",   "6",   "13",   "16"),
  dir   = c("gain", "gain","gain","loss", "loss"),
  lo    = c(123.4,  0,     0,     17.7,   36.8) * 1e6,
  hi    = c(248.9,  93.9,  59.8,  114.4,  90.3) * 1e6
)

MAJ <- fread("results/round_matching_majority.csv",
             colClasses = list(character = "majority_set"))
OLD <- yaml::read_yaml("config/large_clone_comparisons.yaml")

union_len <- function(lo, hi) {
  if (!length(lo)) return(0)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  tot <- 0; cs <- lo[1]; ce <- hi[1]
  for (i in seq_along(lo)[-1]) {
    if (lo[i] > ce) { tot <- tot + (ce - cs); cs <- lo[i]; ce <- hi[i] } else ce <- max(ce, hi[i])
  }
  tot + (ce - cs)
}
split_gt <- function(x) {
  x <- gsub('"', '', x)
  if (is.na(x) || !nzchar(x)) return(character(0))
  sort(trimws(strsplit(x, ",")[[1]]))
}

new <- list(); report <- list()
for (i in seq_len(nrow(MAJ))) {
  s <- MAJ$sample_id[i]; k <- MAJ$best_round[i]
  jpf <- file.path(NB, s, sprintf("joint_post_%d.tsv", k))
  cpf <- file.path(NB, s, sprintf("clone_post_%d.tsv", k))
  if (!file.exists(jpf) || !file.exists(cpf)) next

  jp <- fread(jpf, colClasses = list(character = "CHROM"), showProgress = FALSE,
              select = c("CHROM", "seg", "seg_start", "seg_end", "cnv_state"))
  SEG <- unique(jp)[, .(CHROM = CHROM[1], seg_start = min(seg_start),
                        seg_end = max(seg_end),
                        cnv_state = paste(sort(unique(cnv_state)), collapse = "/")), by = seg]
  setkey(SEG, seg)

  cp <- fread(cpf, showProgress = FALSE, select = c("cell", "clone_opt", "GT_opt"),
              colClasses = list(character = "GT_opt"))
  cl <- cp[, .(n_cells = .N), by = .(clone = clone_opt, GT = GT_opt)]
  gts <- lapply(cl$GT, split_gt); names(gts) <- as.character(cl$clone)
  ncell <- setNames(cl$n_cells, as.character(cl$clone))

  # direct parent of each clone = largest proper subset genotype
  ent <- list()
  for (b in names(gts)) {
    cand <- names(gts)[vapply(names(gts), function(a)
      a != b && all(gts[[a]] %in% gts[[b]]) && length(gts[[a]]) < length(gts[[b]]),
      logical(1))]
    if (!length(cand)) next
    a <- cand[which.max(vapply(cand, function(z) length(gts[[z]]), integer(1)))]
    if (ncell[[a]] < MIN_CELLS || ncell[[b]] < MIN_CELLS) next

    gained <- setdiff(gts[[b]], gts[[a]])
    gs <- SEG[J(gained), nomatch = 0L]
    if (!nrow(gs)) next

    toks <- character(0); segs <- character(0)
    for (e in seq_len(nrow(EV))) {
      pat <- if (EV$dir[e] == "gain") "amp" else "del|loh"
      hit <- gs[CHROM == EV$CHROM[e] & grepl(pat, cnv_state) &
                  seg_end > EV$lo[e] & seg_start < EV$hi[e]]
      if (!nrow(hit)) next
      frac <- union_len(pmax(hit$seg_start, EV$lo[e]), pmin(hit$seg_end, EV$hi[e])) /
        (EV$hi[e] - EV$lo[e])
      if (frac < MIN_ARM) next
      toks <- c(toks, EV$token[e]); segs <- c(segs, hit$seg)
    }
    if (!length(toks)) next
    ent[[sprintf("%s_v_%s_%s", b, a, paste(toks, collapse = "_"))]] <-
      if (length(segs) == 1) segs else as.list(unique(segs))
    report[[length(report) + 1L]] <- data.table(
      sample_id = s, round = k, transition = sprintf("%s_v_%s", b, a),
      events = paste(toks, collapse = ","), segs = paste(unique(segs), collapse = ","),
      cells_a = ncell[[a]], cells_b = ncell[[b]])
  }
  if (length(ent)) new[[s]] <- ent[order(names(ent))]
}

# preserve held-out samples and branch keys verbatim
keep <- names(OLD)[grepl("_branch_", names(OLD)) | !sub("_branch_.*", "", names(OLD)) %in% MAJ$sample_id]
for (kk in keep) new[[kk]] <- OLD[[kk]]

new <- new[order(names(new))]
yaml::write_yaml(new, OUT)

rep <- rbindlist(report)
cat("wrote ", OUT, "\n\n", sep = "")
cat("regenerated entries:", nrow(rep), "across", length(unique(rep$sample_id)), "samples\n")
cat("preserved verbatim  :", length(keep), "keys ->", paste(keep, collapse = " "), "\n\n")
cat("=== per-SCNA sample coverage (regenerated only) ===\n")
for (tk in EV$token) {
  ss <- sort(unique(rep[grepl(tk, events, fixed = TRUE), sample_id]))
  cat(sprintf("%-5s n=%2d : %s\n", tk, length(ss), paste(ss, collapse = " ")))
}
cat("\n=== entries ===\n")
print(rep[order(sample_id)], nrows = 200)
