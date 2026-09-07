#!/usr/bin/env Rscript
# Does genome DOSAGE explain the depth-clone association? (#43 follow-up, part 2)
#
# src/diag_depth_clone_association.R established two things:
#   (a) spearman(nCount_gene, p_opt) > 0 in 36/36 samples, median +0.319 --
#       numbat is less confident about shallow cells, everywhere. Solid.
#   (b) tumor clones separate on depth in 30/30 testable samples, but weakly
#       (median eps2 0.047).
#
# (b) is DIRECTIONALLY UNRESOLVED, and the premise I used to propose it was
# wrong. I claimed sibling subclones "have no biological reason to differ in
# depth". They do: subclones are DEFINED by gains and losses, and genome dosage
# genuinely changes RNA content. A gain-bearing subclone should have more UMIs.
# That mechanism produces exactly the observed signal and is not an artifact.
#
# THIS script resolves the direction. Each clone's GT_opt names the seg_cons
# segments it carries; segs_consensus gives each one a cnv_state, a length and a
# gene content. So every clone has a predictable net dosage, and dosage makes a
# QUANTITATIVE prediction about depth, not just a sign:
#
#   predicted log2 depth fold = log2(1 + sum_s w_s * delta_s)
#     delta: amp +0.5  bamp +1.0  del -0.5  bdel -1.0  loh 0  (copies/2 - 1)
#     w_s  : share of the assessed genome in segment s, by GENE COUNT (primary,
#            since UMIs come from genes) and by BP LENGTH (secondary). Neither is
#            expression-weighted -- the numbat object carries no baseline
#            expression per gene -- so both are approximations to the share of
#            total RNA that segment contributes. They bracket it.
#
# Three outcomes, three different conclusions:
#   * rho(dosage, depth) > 0 AND slope(obs ~ pred) ~ 1  -> (b) is dosage biology
#   * rho > 0 but slope >> 1  -> dosage is real but far too small to account for
#                                the observed spread; something else, i.e. depth,
#                                is driving the partition
#   * rho ~ 0 or < 0          -> dosage does not explain it at all
#
# The delta model is not assumed on faith: segs_consensus$phi_mle is numbat's
# own MEASURED expression fold per segment, and section 0 checks it against the
# assumed values.
#
# The magnitude comparison is reported TWICE: once referenced to the normal
# clone, and once as the spread among TUMOR CLONES ONLY. The first inherits the
# large normal-vs-tumor depth confound (part 1: median log2 ratio +0.75) and so
# overstates the gap; the tumor-only spread is the honest number.
#
# Read-only. Opens the same numbat objects and the same QC database as part 1;
# writes only results/. No numbat rerun, no targets store access.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(numbatHelpers)
})

OUT_CLONE  <- "results/clone_dosage_depth_clones.csv"
OUT_SAMPLE <- "results/clone_dosage_depth_samples.csv"
OUT_PHI    <- "results/clone_dosage_phi_check.csv"
OUT_PDF    <- "results/clone_dosage_depth.pdf"

# relative expression multiplier - 1, i.e. (copies/2) - 1
DELTA <- c(neu = 0, amp = 0.5, bamp = 1.0, del = -0.5, bdel = -1.0, loh = 0)

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2)
    return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

rng <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else diff(range(x)) }

rds <- sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) rds <- rds[grepl(paste(args, collapse = "|"), rds)]
cat("samples:", length(rds), "\n")

clone_rows <- list(); samp_rows <- list(); phi_rows <- list()

for (f in rds) {
  s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) { cat("  skip", s, "(unreadable)\n"); next }

  cp <- tryCatch(as.data.frame(nb$clone_post), error = function(e) NULL)
  sc <- tryCatch(as.data.frame(nb$segs_consensus), error = function(e) NULL)
  gtf <- tryCatch(as.data.frame(nb$gtf), error = function(e) NULL)
  rm(nb); invisible(gc(FALSE))  # numbat objects are GB-scale; do not let them stack
  if (is.null(cp) || !nrow(cp) || is.null(sc) || !nrow(sc) ||
      !all(c("seg_cons", "seg_length", "seg_start", "seg_end") %in% names(sc))) {
    cat("  skip", s, "(no clone_post / segs_consensus)\n"); next
  }
  if (!"GT_opt" %in% names(cp)) { cat("  skip", s, "(no GT_opt)\n"); next }

  q <- tryCatch(get_sample_cell_depth(s), error = function(e) NULL)
  if (is.null(q) || !nrow(q)) { cat("  skip", s, "(no cell_qc_values)\n"); next }

  st <- if ("cnv_state_post" %in% names(sc)) "cnv_state_post" else "cnv_state"
  sc$chr_i <- suppressWarnings(as.integer(as.character(sc$CHROM)))

  # -- weights -----------------------------------------------------------------
  # Denominator = the assessed genome. segs_consensus lists a physical segment
  # ONCE PER ALTERNATIVE STATE (1c/amp and 1d/del are the same coordinates), so
  # dedup on coordinates or the total is inflated by every branching segment.
  phys <- sc |> distinct(chr_i, seg_start, seg_end, .keep_all = TRUE)
  total_len <- sum(phys$seg_length, na.rm = TRUE)
  if (!is.finite(total_len) || total_len <= 0) { cat("  skip", s, "(bad lengths)\n"); next }

  # gene share: count gtf genes by midpoint inside each segment. gtf CHROM is
  # integer, segs_consensus CHROM is a factor -- coerce or nothing ever matches.
  seg_wg <- NULL
  if (!is.null(gtf) && all(c("gene_start", "gene_end", "CHROM") %in% names(gtf))) {
    gmid <- (gtf$gene_start + gtf$gene_end) / 2
    gchr <- suppressWarnings(as.integer(as.character(gtf$CHROM)))
    ngene <- vapply(seq_len(nrow(sc)), function(i)
      sum(gchr == sc$chr_i[i] & gmid >= sc$seg_start[i] & gmid <= sc$seg_end[i],
          na.rm = TRUE), integer(1))
    nphys <- vapply(seq_len(nrow(phys)), function(i)
      sum(gchr == phys$chr_i[i] & gmid >= phys$seg_start[i] & gmid <= phys$seg_end[i],
          na.rm = TRUE), integer(1))
    tot_gene <- sum(nphys)
    if (tot_gene > 0) seg_wg <- setNames(ngene / tot_gene, sc$seg_cons)
  }
  if (is.null(seg_wg)) { cat("  note", s, "(no gtf; gene weight unavailable)\n") }

  seg_wl <- setNames(sc$seg_length / total_len, sc$seg_cons)
  seg_d  <- setNames(unname(DELTA[sc[[st]]]), sc$seg_cons)

  if ("phi_mle" %in% names(sc))
    phi_rows[[s]] <- tibble::tibble(sample_id = s, seg_cons = sc$seg_cons,
                                    cnv_state = sc[[st]], phi_mle = sc$phi_mle,
                                    assumed_fold = 1 + unname(DELTA[sc[[st]]]))

  gt <- unique(cp[, c("clone_opt", "GT_opt"), drop = FALSE])
  gt <- gt[!duplicated(gt$clone_opt), ]

  tok_tot <- 0L; tok_hit <- 0L
  dosage_with <- function(w) vapply(seq_len(nrow(gt)), function(i) {
    g <- gt$GT_opt[i]
    if (is.na(g) || !nzchar(trimws(g))) return(0)          # normal clone
    tk <- trimws(strsplit(g, ",")[[1]]); tk <- tk[nzchar(tk)]
    hit <- tk %in% names(w)
    if (!any(hit)) return(NA_real_)
    sum(w[tk[hit]] * seg_d[tk[hit]], na.rm = TRUE)
  }, numeric(1))

  for (i in seq_len(nrow(gt))) {                            # token QC, once
    g <- gt$GT_opt[i]
    if (is.na(g) || !nzchar(trimws(g))) next
    tk <- trimws(strsplit(g, ",")[[1]]); tk <- tk[nzchar(tk)]
    tok_tot <- tok_tot + length(tk); tok_hit <- tok_hit + sum(tk %in% names(seg_wl))
  }

  gt$net_dosage_len  <- dosage_with(seg_wl)
  gt$net_dosage_gene <- if (!is.null(seg_wg)) dosage_with(seg_wg) else NA_real_
  # primary = gene-weighted where available, else length
  gt$net_dosage <- ifelse(is.finite(gt$net_dosage_gene),
                          gt$net_dosage_gene, gt$net_dosage_len)
  gt$is_normal  <- !is.na(gt$GT_opt) & trimws(gt$GT_opt) == ""

  m <- cp[, base::intersect(c("cell", "clone_opt", "p_opt"), names(cp)), drop = FALSE] |>
    inner_join(as.data.frame(q), by = "cell") |>
    inner_join(gt[, c("clone_opt", "net_dosage", "net_dosage_len",
                      "net_dosage_gene", "is_normal")], by = "clone_opt")
  if (nrow(m) < 50) { cat("  skip", s, "(", nrow(m), "joined cells )\n"); next }

  cl <- m |> group_by(clone_opt, net_dosage, net_dosage_len, net_dosage_gene, is_normal) |>
    summarise(n_cells = n(), median_depth = median(nCount_gene), .groups = "drop") |>
    mutate(sample_id = s)

  # Reference clone for the observed fold: the normal clone if numbat found one,
  # otherwise the tumor clone closest to neutral dosage.
  ref <- if (any(cl$is_normal)) cl$median_depth[cl$is_normal][1]
         else cl$median_depth[which.min(abs(cl$net_dosage))]
  cl <- cl |> mutate(
    obs_log2fc  = log2(median_depth / ref),
    pred_log2fc = ifelse(is.finite(net_dosage) & net_dosage > -1,
                         log2(1 + net_dosage), NA_real_))
  clone_rows[[s]] <- cl

  tu <- cl[!cl$is_normal, ]
  samp_rows[[s]] <- tibble::tibble(
    sample_id       = s,
    n_cells         = nrow(m),
    n_clones        = nrow(cl),
    n_tumor_clones  = nrow(tu),
    gt_token_match  = if (tok_tot) tok_hit / tok_tot else NA_real_,
    dosage_range    = rng(cl$net_dosage),
    # clone-level: does depth track dosage across clones?
    rho_clone_all   = sp(cl$net_dosage, cl$median_depth),
    rho_clone_tumor = sp(tu$net_dosage, tu$median_depth),
    # cell-level: same question with all the cells behind it
    rho_cell_all      = sp(m$net_dosage, m$nCount_gene),
    rho_cell_tumor    = sp(m$net_dosage[!m$is_normal], m$nCount_gene[!m$is_normal]),
    rho_cell_tumor_len = sp(m$net_dosage_len[!m$is_normal], m$nCount_gene[!m$is_normal]),
    # magnitude, normal-referenced (inherits the normal/tumor confound)
    max_abs_pred    = suppressWarnings(max(abs(tu$pred_log2fc), na.rm = TRUE)),
    max_abs_obs     = suppressWarnings(max(abs(tu$obs_log2fc),  na.rm = TRUE)),
    # magnitude, TUMOR CLONES ONLY -- the honest comparison
    tumor_pred_range = rng(tu$pred_log2fc),
    tumor_obs_range  = rng(tu$obs_log2fc))

  cat(sprintf("  %s  cells %5d  clones %d (tumor %d)  GT %.0f%%  dose[gene] %+.4f..%+.4f  rho_cell_tumor %s\n",
              s, nrow(m), nrow(cl), nrow(tu), 100 * (if (tok_tot) tok_hit/tok_tot else NA),
              min(cl$net_dosage, na.rm = TRUE), max(cl$net_dosage, na.rm = TRUE),
              ifelse(is.na(samp_rows[[s]]$rho_cell_tumor), "NA",
                     sprintf("%+.3f", samp_rows[[s]]$rho_cell_tumor))))
}

stopifnot(length(samp_rows) > 0)
S <- bind_rows(samp_rows)
C <- bind_rows(clone_rows) |> relocate(sample_id)
P <- bind_rows(phi_rows)
S$max_abs_pred[!is.finite(S$max_abs_pred)] <- NA_real_
S$max_abs_obs[!is.finite(S$max_abs_obs)]   <- NA_real_

readr::write_csv(S, OUT_SAMPLE); readr::write_csv(C, OUT_CLONE)
if (nrow(P)) readr::write_csv(P, OUT_PHI)

sgn <- function(x, lab) {
  x <- x[is.finite(x)]
  if (!length(x)) { cat(sprintf("  %-22s  (none testable)\n", lab)); return(invisible()) }
  cat(sprintf("  %-22s median %+.3f   positive in %2d/%2d samples   (sign test p = %.3g)\n",
              lab, median(x), sum(x > 0), length(x),
              binom.test(sum(x > 0), length(x))$p.value))
}

cat("\n================ RESULT ================\n")
cat(sprintf("samples %d   clones %d   GT_opt token match rate: median %.1f%%, min %.1f%%\n",
            nrow(S), nrow(C), 100*median(S$gt_token_match, na.rm=TRUE),
            100*min(S$gt_token_match, na.rm=TRUE)))

cat("\n-- 0. is the dosage model right? numbat's own measured phi_mle by state --\n")
if (nrow(P)) {
  pk <- P |> filter(is.finite(phi_mle)) |> group_by(cnv_state) |>
    summarise(n = n(), assumed = first(assumed_fold),
              median_phi = median(phi_mle), q25 = quantile(phi_mle,.25),
              q75 = quantile(phi_mle,.75), .groups = "drop")
  print(as.data.frame(pk), row.names = FALSE, digits = 3)
}

cat("\n-- 1. does depth track dosage? (positive = gain-bearing clones deeper) --\n")
sgn(S$rho_clone_all,      "clone, all")
sgn(S$rho_clone_tumor,    "clone, tumor only")
sgn(S$rho_cell_all,       "cell,  all")
sgn(S$rho_cell_tumor,     "cell,  tumor only")
sgn(S$rho_cell_tumor_len, "cell,  tumor (bp wt)")

cat("\n-- 2. magnitude: observed depth spread vs what dosage predicts --\n")
Ct <- C |> filter(!is_normal, is.finite(pred_log2fc), is.finite(obs_log2fc))
fit <- if (nrow(Ct) >= 10) lm(obs_log2fc ~ pred_log2fc, data = Ct) else NULL
if (!is.null(fit)) {
  co <- summary(fit)$coefficients
  cat(sprintf("  lm(obs ~ pred), %d tumor clones, normal-referenced: slope %+.2f (SE %.2f, p %.3g), R2 %.3f\n",
              nrow(Ct), co[2,1], co[2,2], co[2,4], summary(fit)$r.squared))
  cat("    slope ~ 1 => dosage accounts for the depth differences\n")
  cat("    slope >> 1 => observed differences far exceed dosage; depth is driving\n")
}
mr <- S |> filter(is.finite(max_abs_pred), is.finite(max_abs_obs), max_abs_pred > 0)
if (nrow(mr))
  cat(sprintf("  normal-referenced max|obs|/max|pred|: median %.0fx  (obs %.3f vs pred %.4f)  [n=%d]\n",
              median(mr$max_abs_obs / mr$max_abs_pred),
              median(mr$max_abs_obs), median(mr$max_abs_pred), nrow(mr)))
tr <- S |> filter(is.finite(tumor_pred_range), is.finite(tumor_obs_range), tumor_pred_range > 0)
if (nrow(tr)) {
  cat(sprintf("  TUMOR-ONLY spread  obs/pred:        median %.0fx  (obs %.3f vs pred %.4f)  [n=%d]\n",
              median(tr$tumor_obs_range / tr$tumor_pred_range),
              median(tr$tumor_obs_range), median(tr$tumor_pred_range), nrow(tr)))
  cat(sprintf("    obs spread exceeds pred in %d/%d samples\n",
              sum(tr$tumor_obs_range > tr$tumor_pred_range), nrow(tr)))
}

cat("\n-- 3. how much dosage variation is there to detect? --\n")
cat(sprintf("  net dosage range within sample: median %.4f, max %.4f\n",
            median(S$dosage_range, na.rm=TRUE), max(S$dosage_range, na.rm=TRUE)))
cat(sprintf("  = a median predicted depth difference of %.2f%% between extreme clones\n",
            100*median(S$dosage_range, na.rm=TRUE)))

cat("\n-- per-sample --\n")
print(as.data.frame(S |> arrange(desc(abs(rho_cell_tumor))) |>
  select(sample_id, n_cells, n_tumor_clones, dosage_range, rho_cell_tumor,
         rho_clone_tumor, tumor_pred_range, tumor_obs_range)),
  row.names = FALSE, digits = 3)

# ---- figure -----------------------------------------------------------------
p1 <- ggplot(Ct, aes(pred_log2fc, obs_log2fc)) +
  geom_abline(slope = 1, intercept = 0, colour = "#2980b9", linewidth = .7) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .3) +
  geom_point(aes(size = n_cells), alpha = .6, colour = "#c0392b") +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = .6) +
  scale_size_area(max_size = 5) +
  labs(title = "Observed clone depth difference vs what genome dosage predicts",
       subtitle = paste("blue = dosage fully explains it (slope 1); black = fit.",
                        "Tumor clones, referenced to the normal clone.",
                        "\nNote the x scale: predicted dosage effects are ~1%, observed differences are tens of percent."),
       x = "predicted log2 fold from net dosage (gene-weighted)",
       y = "observed log2 fold in nCount_gene") +
  theme_bw(base_size = 9)

p2 <- ggplot(S, aes(reorder(sample_id, rho_cell_tumor), rho_cell_tumor)) +
  geom_col(fill = "#8e44ad") + coord_flip() +
  geom_hline(yintercept = 0, colour = "grey30") +
  labs(title = "Within-sample correlation of clone dosage with cell depth",
       subtitle = "tumor cells only; positive = gain-bearing clones are deeper, as dosage predicts",
       x = NULL, y = "spearman(net dosage, nCount_gene)") +
  theme_bw(base_size = 8)

p3 <- ggplot(tr, aes(tumor_pred_range, tumor_obs_range)) +
  geom_abline(slope = 1, intercept = 0, colour = "#2980b9", linewidth = .7) +
  geom_point(colour = "#c0392b", size = 2, alpha = .7) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Tumor-clone depth spread vs dosage-predicted spread",
       subtitle = "normal clone excluded, so the normal-vs-tumor confound cannot contribute. Points above the line = depth varies more than dosage can explain.",
       x = "predicted log2 spread across tumor clones",
       y = "observed log2 spread across tumor clones") +
  theme_bw(base_size = 9)

p4 <- ggplot(C, aes(net_dosage, median_depth, colour = is_normal)) +
  geom_point(size = 1) +
  scale_colour_manual(values = c(`TRUE` = "grey60", `FALSE` = "#c0392b"),
                      name = "normal clone") +
  scale_y_log10() +
  facet_wrap(~ sample_id, scales = "free", ncol = 6) +
  labs(title = "Clone median depth vs net genome dosage",
       x = "net dosage (gene-share weighted)", y = "median nCount_gene") +
  theme_bw(base_size = 6) + theme(legend.position = "top")

grDevices::pdf(OUT_PDF, width = 12, height = 9)
print(p1); print(p2); print(p3); print(p4)
invisible(grDevices::dev.off())

cat("\nwrote", OUT_SAMPLE, "/", OUT_CLONE, "/", OUT_PHI, "/", OUT_PDF, "\nDIAG DONE\n")
