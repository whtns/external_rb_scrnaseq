# Two-clone comparisons in which a LARGE 2p or 6p gain is the dominant change.
#
# The use case is attributing expression differences to a single SCNA. That
# needs a pair of clones from the same tumour differing by the 2p/6p gain and
# by as little else as possible -- a pair that also differs by 1q, 16q and 13q
# tells you nothing about 2p specifically.
#
# Three requirements, all applied here:
#   1. CLEAN     few other segments differ between the two clones
#   2. LARGE     the 2p/6p gain itself spans a substantial part of the arm
#   3. MYCN      2p gains must overlap MYCN (chr2:15,940,550-15,947,007, hg38);
#                a 2p gain that misses MYCN is not the event of interest
#
# LABEL SPACE WARNING. Clone genotypes (`GT_opt` in clone_post) use the segment
# labels of `joint_post`, NOT those of `segs_consensus` -- in SRX10031193 round
# 2 the GT segments 2f/2i/2m/6c exist only in joint_post, while segs_consensus
# calls the same region 2e/2g/2k. Joining GT labels to segs_consensus silently
# drops most of the genotype. Coordinates here come from joint_post.
#
# Reads only. Nothing writes into output/numbat_sridhar/.
#
# Outputs:
#   results/clone_pairs_2p_6p_all.csv    every candidate pair, ranked
#   results/clone_pairs_2p_6p_clean.csv  the shortlist meeting all three
#
# Usage: Rscript src/find_clean_2p_6p_clone_pairs.R

suppressPackageStartupMessages(library(data.table))

NB_DIR    <- "output/numbat_sridhar"
MYCN      <- c(15940550, 15947007)
MIN_CELLS <- 20      # per clone, below which a diffex comparison is not viable
MIN_MB    <- 20      # "large": focal gain extent inside the arm window

# Arm windows, as in numbat_helpers::numbat_rb_arm_frac (hg38 centromeres).
WIN <- data.table(event = c("2p_gain", "6p_gain"), CHROM = c("2", "6"),
                  lo = c(0, 0), hi = c(93.9e6, 59.8e6))

MAJ <- fread("results/round_matching_majority.csv",
             colClasses = list(character = "majority_set"))

union_len <- function(lo, hi) {
  if (!length(lo)) return(0)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  tot <- 0; cs <- lo[1]; ce <- hi[1]
  for (i in seq_along(lo)[-1]) {
    if (lo[i] > ce) { tot <- tot + (ce - cs); cs <- lo[i]; ce <- hi[i] }
    else ce <- max(ce, hi[i])
  }
  tot + (ce - cs)
}

split_gt <- function(x) {
  x <- gsub('"', '', x)
  if (is.na(x) || !nzchar(x)) return(character(0))
  sort(trimws(strsplit(x, ",")[[1]]))
}

rows <- list()
for (i in seq_len(nrow(MAJ))) {
  s <- MAJ$sample_id[i]; k <- MAJ$best_round[i]

  jpf <- file.path(NB_DIR, s, sprintf("joint_post_%d.tsv", k))
  cpf <- file.path(NB_DIR, s, sprintf("clone_post_%d.tsv", k))
  if (!file.exists(jpf) || !file.exists(cpf)) { warning("missing files for ", s); next }

  jp <- fread(jpf, colClasses = list(character = "CHROM"), showProgress = FALSE,
              select = c("CHROM", "seg", "seg_start", "seg_end", "cnv_state"))
  SEG <- unique(jp)
  # one row per label; a label carrying >1 state would break the join
  SEG <- SEG[, .(CHROM = CHROM[1], seg_start = min(seg_start), seg_end = max(seg_end),
                 cnv_state = paste(sort(unique(cnv_state)), collapse = "/")), by = seg]
  setkey(SEG, seg)

  cp <- fread(cpf, showProgress = FALSE, select = c("cell", "clone_opt", "GT_opt"),
              colClasses = list(character = "GT_opt"))
  cl <- cp[, .(n_cells = .N), by = .(clone = clone_opt, GT = GT_opt)]
  if (nrow(cl) < 2) next
  gts <- lapply(cl$GT, split_gt)
  names(gts) <- as.character(cl$clone)

  for (a in seq_len(nrow(cl))) for (b in seq_len(nrow(cl))) {
    if (a == b) next
    A <- gts[[a]]; B <- gts[[b]]
    gained <- setdiff(B, A); lost <- setdiff(A, B)
    if (!length(gained)) next
    nested <- length(lost) == 0

    gs <- SEG[J(gained), nomatch = 0L]
    if (!nrow(gs)) next

    for (w in seq_len(nrow(WIN))) {
      ev <- WIN$event[w]
      foc <- gs[CHROM == WIN$CHROM[w] & grepl("amp", cnv_state) &
                  seg_end > WIN$lo[w] & seg_start < WIN$hi[w]]
      if (!nrow(foc)) next

      flo <- pmax(foc$seg_start, WIN$lo[w]); fhi <- pmin(foc$seg_end, WIN$hi[w])
      focal_mb <- union_len(flo, fhi) / 1e6
      # MYCN is a chr2 locus, so the test is only meaningful for 2p. Left as a
      # plain overlap test it would return TRUE for any chr6 segment spanning
      # 15.9 Mb -- a coordinate coincidence, not a fact about the sample.
      hits_mycn <- if (WIN$CHROM[w] == "2")
        any(foc$seg_end > MYCN[1] & foc$seg_start < MYCN[2]) else NA

      # everything else that differs between the two clones
      oth <- rbind(gs[!seg %in% foc$seg], SEG[J(lost), nomatch = 0L])
      rows[[length(rows) + 1L]] <- data.table(
        sample_id = s, round = k, event = ev,
        clone_a = cl$clone[a], clone_b = cl$clone[b],
        n_cells_a = cl$n_cells[a], n_cells_b = cl$n_cells[b],
        relation      = if (nested) "nested" else "sibling",
        focal_segs    = paste(foc$seg, collapse = ","),
        focal_mb      = round(focal_mb, 1),
        focal_arm_frac = round(focal_mb * 1e6 / (WIN$hi[w] - WIN$lo[w]), 3),
        focal_state   = paste(sort(unique(foc$cnv_state)), collapse = ","),
        # Unclipped extent. A "2p gain" whose segment runs to 170 Mb is a
        # whole-chromosome-2 gain: 2q rides along and cannot be separated from
        # 2p/MYCN by this comparison. Worth knowing before using the pair.
        focal_full_mb = round((max(foc$seg_end) - min(foc$seg_start)) / 1e6, 1),
        beyond_arm_mb = round(max(0, (max(foc$seg_end) - WIN$hi[w])) / 1e6, 1),
        whole_chrom   = max(foc$seg_end) > WIN$hi[w] * 1.25,
        mycn          = hits_mycn,
        n_other_seg   = nrow(oth),
        n_other_chrom = length(unique(oth$CHROM)),
        other_mb      = round(sum(oth$seg_end - oth$seg_start) / 1e6, 1),
        other_detail  = if (nrow(oth))
          paste(sprintf("%s:%s(%s)", oth$CHROM, oth$seg, oth$cnv_state), collapse = ";") else ""
      )
    }
  }
}

res <- rbindlist(rows)
res[, viable := n_cells_a >= MIN_CELLS & n_cells_b >= MIN_CELLS]
res[, large  := focal_mb >= MIN_MB]
# 2p must hit MYCN; 6p has no such requirement
res[, passes := viable & large & (event == "6p_gain" | mycn)]
setorder(res, event, n_other_chrom, n_other_seg, -focal_mb)
fwrite(res, "results/clone_pairs_2p_6p_all.csv")

clean <- res[passes == TRUE]
fwrite(clean, "results/clone_pairs_2p_6p_clean.csv")

show <- c("sample_id", "event", "clone_a", "clone_b", "n_cells_a", "n_cells_b",
          "relation", "focal_mb", "focal_arm_frac", "mycn",
          "n_other_seg", "n_other_chrom", "other_mb")
cat("\n===== candidates passing viability + size (+ MYCN for 2p) =====\n")
cat("ranked by fewest other changes; n_other_seg == 0 means the ONLY difference\n")
for (ev in c("2p_gain", "6p_gain")) {
  cat("\n--- ", ev, " ---\n", sep = "")
  d <- clean[event == ev]
  if (!nrow(d)) { cat("(none)\n"); next }
  print(head(d[, ..show], 25))
}
# One row per (sample, event): the cleanest viable pair available for it.
best <- clean[order(event, sample_id, n_other_chrom, n_other_seg, -focal_mb)][
  , .SD[1], by = .(sample_id, event)]
setorder(best, event, n_other_chrom, n_other_seg, -focal_mb)
fwrite(best, "results/clone_pairs_2p_6p_best_per_sample.csv")
cat("\n===== BEST PAIR PER SAMPLE =====\n")
for (ev in c("2p_gain", "6p_gain")) {
  cat("\n--- ", ev, " ---\n", sep = "")
  d <- best[event == ev]
  if (!nrow(d)) { cat("(none)\n"); next }
  print(d[, .(sample_id, clone_a, clone_b, n_cells_a, n_cells_b, relation,
              focal_mb, focal_arm_frac, whole_chrom, mycn, n_other_seg,
              n_other_chrom, other_mb, other_detail)])
}

cat("\n===== why candidates were excluded =====\n")
print(res[, .(n = .N,
              fail_cells = sum(!viable),
              fail_size  = sum(viable & !large),
              fail_mycn  = sum(viable & large & event == "2p_gain" & !mycn),
              pass       = sum(passes)), by = event])
