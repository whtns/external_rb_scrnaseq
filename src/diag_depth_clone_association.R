#!/usr/bin/env Rscript
# Are numbat clones an artifact of read depth? (github #43 follow-up)
#
# Everything #43 delivered is BETWEEN-sample and descriptive: a marginal depth
# panel per summary, plus one weak sample-level correlation
# (spearman(median depth, n_rb_events) = 0.282, p = 0.10, n = 35; it falls to
# 0.219 / p = 0.21 once the single broken sample SRX10031191 is dropped).
#
# "Are clones artifacts of depth" is a WITHIN-sample, PER-CELL question, and it
# has never been asked. The data to ask it already exists and joins cleanly:
# nb$clone_post carries (cell, clone_opt, p_opt) and cell_qc_values carries
# (cell, nCount_gene, ...) on the same barcodes -- 4085/4085 on the pilot.
#
# Pilot (SRX10264519): Kruskal-Wallis depth ~ clone p = 2.5e-10, and
# spearman(depth, p_opt) = +0.537. Both are what a depth artifact would look
# like. Both are ALSO what real biology looks like, which is the point of the
# controls below.
#
# THE CONFOUND, and how this handles it. numbat's clone 1 has an empty GT_opt --
# it is the normal/diploid clone, verified on SRX10264519 and SRX22868104.
# Normal and malignant cells differ in RNA content for reasons that have nothing
# to do with sequencing artifact, so a depth split across ALL clones is expected
# and uninformative. The discriminating test is the same statistic computed over
# the TUMOR clones only: sibling subclones of one tumor have no biological
# reason to differ systematically in depth, so if they still separate, depth is
# driving the partition rather than the reverse.
#
# Reports, per sample:
#   kw_p_all      Kruskal-Wallis depth ~ clone, every clone      (expected sig.)
#   kw_p_tumor    the same over tumor clones only                (the real test)
#   eps2_*        epsilon-squared effect size for each, so a significant p on
#                 30k cells is not mistaken for a large effect
#   rho_p_opt     spearman(depth, p_opt) -- is numbat less SURE on shallow cells
#   rho_p_cnv     spearman(depth, mean per-cell p_cnv over called segments)
#
# BH-adjusted across samples. Read-only: opens the numbat objects and the QC
# database, writes only results/ and a figure.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(numbatHelpers)
})

OUT_CSV   <- "results/depth_clone_association.csv"
OUT_CLONE <- "results/depth_clone_per_clone.csv"
OUT_CELL  <- "results/depth_clone_cells.csv.gz"
OUT_PDF   <- "results/depth_clone_association.pdf"

eps2 <- function(H, k, n) if (is.na(H) || n <= k) NA_real_ else (H - k + 1) / (n - k)

kw <- function(depth, grp) {
  grp <- factor(grp)
  if (nlevels(grp) < 2) return(list(p = NA_real_, eps2 = NA_real_, k = nlevels(grp)))
  t <- suppressWarnings(kruskal.test(depth ~ grp))
  list(p = unname(t$p.value),
       eps2 = eps2(unname(t$statistic), nlevels(grp), length(depth)),
       k = nlevels(grp))
}

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20 || length(unique(y[ok])) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

rds <- sort(Sys.glob("output/numbat_sridhar/SRX*_numbat.rds"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) rds <- rds[grepl(paste(args, collapse = "|"), rds)]
cat("samples:", length(rds), "\n")

cells <- list(); per_clone <- list(); per_sample <- list()

for (f in rds) {
  s  <- stringr::str_extract(basename(f), "SRX[0-9]+")
  nb <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(nb)) { cat("  skip", s, "(unreadable)\n"); next }

  cp <- tryCatch(nb$clone_post, error = function(e) NULL)
  if (is.null(cp) || !nrow(cp)) { cat("  skip", s, "(no clone_post)\n"); next }

  q <- tryCatch(get_sample_cell_depth(s), error = function(e) NULL)
  if (is.null(q) || !nrow(q)) { cat("  skip", s, "(no cell_qc_values)\n"); next }

  # per-cell mean posterior over the segments numbat actually called
  jp <- tryCatch(nb$joint_post, error = function(e) NULL)
  pc <- NULL
  if (!is.null(jp) && all(c("cell", "seg", "p_cnv") %in% names(jp))) {
    sc <- tryCatch(nb$segs_consensus, error = function(e) NULL)
    if (!is.null(sc)) {
      st <- if ("cnv_state_post" %in% names(sc)) "cnv_state_post" else "cnv_state"
      called <- sc$seg[sc[[st]] != "neu"]
      if (length(called))
        pc <- as.data.frame(jp) |> filter(seg %in% called) |> group_by(cell) |>
          summarise(mean_p_cnv = mean(p_cnv, na.rm = TRUE), .groups = "drop")
    }
  }

  # clone_post is a data.table; `[` with a character j and drop= does not
  # subset columns the way a data.frame does (it returned a bare character
  # vector, which inner_join then rejected). Coerce first, and use
  # base::intersect -- numbatHelpers masks intersect with GenomeInfoDb's.
  cp   <- as.data.frame(cp)
  keep <- base::intersect(c("cell", "clone_opt", "p_opt", "GT_opt"), names(cp))
  m <- cp[, keep, drop = FALSE] |>
    inner_join(as.data.frame(q), by = "cell")
  if (!is.null(pc)) m <- left_join(m, pc, by = "cell")
  if (!"mean_p_cnv" %in% names(m)) m$mean_p_cnv <- NA_real_
  if (!"GT_opt" %in% names(m))     m$GT_opt     <- NA_character_
  if (nrow(m) < 50) { cat("  skip", s, "(", nrow(m), "joined cells )\n"); next }

  # normal clone = empty genotype string; everything else is tumor
  m$is_normal <- !is.na(m$GT_opt) & trimws(m$GT_opt) == ""
  m$sample_id <- s
  cat(sprintf("  %s  cells %5d/%-5d  clones %d  normal-clone cells %d\n",
              s, nrow(m), nrow(cp), dplyr::n_distinct(m$clone_opt), sum(m$is_normal)))

  a <- kw(m$nCount_gene, m$clone_opt)
  tu <- m[!m$is_normal, ]
  t  <- if (nrow(tu) >= 50) kw(tu$nCount_gene, tu$clone_opt)
        else list(p = NA_real_, eps2 = NA_real_, k = NA_integer_)

  per_sample[[s]] <- tibble::tibble(
    sample_id      = s,
    n_cells        = nrow(m),
    n_clones       = a$k,
    n_tumor_clones = t$k,
    n_normal_cells = sum(m$is_normal),
    median_depth   = median(m$nCount_gene),
    med_depth_normal = if (any(m$is_normal)) median(m$nCount_gene[m$is_normal]) else NA_real_,
    med_depth_tumor  = if (any(!m$is_normal)) median(m$nCount_gene[!m$is_normal]) else NA_real_,
    kw_p_all       = a$p,  eps2_all   = a$eps2,
    kw_p_tumor     = t$p,  eps2_tumor = t$eps2,
    rho_p_opt      = sp(m$nCount_gene, m$p_opt),
    rho_p_cnv      = sp(m$nCount_gene, m$mean_p_cnv))

  per_clone[[s]] <- m |> group_by(sample_id, clone_opt) |>
    summarise(is_normal = any(is_normal), n_cells = n(),
              median_depth = median(nCount_gene),
              q25 = quantile(nCount_gene, .25), q75 = quantile(nCount_gene, .75),
              median_p_opt = median(p_opt, na.rm = TRUE),
              median_pct_mt = median(percent_mt, na.rm = TRUE), .groups = "drop")

  cells[[s]] <- m |> select(sample_id, cell, clone_opt, is_normal,
                            nCount_gene, nFeature_gene, percent_mt, p_opt, mean_p_cnv)
}

stopifnot(length(per_sample) > 0)
S <- bind_rows(per_sample); C <- bind_rows(per_clone); X <- bind_rows(cells)

S <- S |> mutate(kw_q_all   = p.adjust(kw_p_all,   "BH"),
                 kw_q_tumor = p.adjust(kw_p_tumor, "BH"))

readr::write_csv(S, OUT_CSV); readr::write_csv(C, OUT_CLONE)
readr::write_csv(X, OUT_CELL)

cat("\n================ RESULT ================\n")
cat(sprintf("samples analysed: %d   cells: %d\n", nrow(S), nrow(X)))

cat("\n-- all clones (normal vs tumor split expected; NOT the test) --\n")
cat(sprintf("  BH q < 0.05: %d / %d   median eps2 = %.4f\n",
            sum(S$kw_q_all < 0.05, na.rm = TRUE), sum(!is.na(S$kw_q_all)),
            median(S$eps2_all, na.rm = TRUE)))

cat("\n-- TUMOR CLONES ONLY (the discriminating test) --\n")
tt <- S[!is.na(S$kw_q_tumor), ]
cat(sprintf("  testable samples: %d\n", nrow(tt)))
cat(sprintf("  BH q < 0.05: %d / %d   median eps2 = %.4f\n",
            sum(tt$kw_q_tumor < 0.05), nrow(tt), median(tt$eps2_tumor, na.rm = TRUE)))
cat(sprintf("  eps2 > 0.06 (conventionally 'moderate'): %d\n",
            sum(tt$eps2_tumor > 0.06, na.rm = TRUE)))

cat("\n-- is numbat less certain on shallow cells? --\n")
cat(sprintf("  spearman(depth, p_opt)     median %+.3f   positive in %d/%d samples\n",
            median(S$rho_p_opt, na.rm = TRUE), sum(S$rho_p_opt > 0, na.rm = TRUE),
            sum(!is.na(S$rho_p_opt))))
cat(sprintf("  spearman(depth, mean p_cnv) median %+.3f   positive in %d/%d samples\n",
            median(S$rho_p_cnv, na.rm = TRUE), sum(S$rho_p_cnv > 0, na.rm = TRUE),
            sum(!is.na(S$rho_p_cnv))))

cat("\n-- normal vs tumor clone depth (the expected confound) --\n")
nz <- S[!is.na(S$med_depth_normal) & !is.na(S$med_depth_tumor), ]
cat(sprintf("  samples with both: %d   tumor deeper in %d\n",
            nrow(nz), sum(nz$med_depth_tumor > nz$med_depth_normal)))
if (nrow(nz)) cat(sprintf("  median(log2 tumor/normal depth ratio) = %+.3f\n",
                          median(log2(nz$med_depth_tumor / nz$med_depth_normal))))

cat("\n-- strongest tumor-only effects --\n")
print(as.data.frame(tt |> arrange(desc(eps2_tumor)) |>
  select(sample_id, n_cells, n_tumor_clones, eps2_tumor, kw_q_tumor, rho_p_opt) |>
  head(10)), row.names = FALSE)

# ---- figure -----------------------------------------------------------------
lab <- S |> mutate(
  ttl = sprintf("%s  (eps2_tumor %s)", sample_id,
                ifelse(is.na(eps2_tumor), "NA", sprintf("%.3f", eps2_tumor)))) |>
  select(sample_id, ttl)

X2 <- X |> left_join(lab, by = "sample_id") |>
  mutate(clone = factor(clone_opt),
         kind  = ifelse(is_normal, "normal", "tumor"))

p1 <- ggplot(X2, aes(clone, nCount_gene, fill = kind)) +
  geom_boxplot(outlier.size = .2, linewidth = .3) +
  scale_y_log10() +
  scale_fill_manual(values = c(normal = "grey70", tumor = "#c0392b")) +
  facet_wrap(~ ttl, scales = "free", ncol = 5) +
  labs(title = "Cell depth by numbat clone",
       subtitle = paste("normal clone = empty GT_opt. The test is whether the TUMOR clones",
                        "separate on depth;\na normal-vs-tumor difference is expected and not evidence of artifact."),
       x = "clone_opt", y = "nCount_gene (log10)") +
  theme_bw(base_size = 7) + theme(legend.position = "top")

p2 <- ggplot(tt, aes(reorder(sample_id, eps2_tumor), eps2_tumor,
                     fill = kw_q_tumor < 0.05)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 0.06, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey75"),
                    name = "BH q < 0.05") +
  labs(title = "Depth ~ clone effect size, TUMOR CLONES ONLY",
       subtitle = "epsilon-squared from Kruskal-Wallis; dashed line = conventional 'moderate' (0.06)",
       x = NULL, y = "epsilon-squared") +
  theme_bw(base_size = 8)

p3 <- ggplot(S, aes(reorder(sample_id, rho_p_opt), rho_p_opt)) +
  geom_col(fill = "#2c3e50") + coord_flip() +
  geom_hline(yintercept = 0, colour = "grey40") +
  labs(title = "Is numbat less certain about shallow cells?",
       subtitle = "spearman(nCount_gene, p_opt); positive = shallow cells get lower-confidence clone calls",
       x = NULL, y = "rho") +
  theme_bw(base_size = 8)

grDevices::pdf(OUT_PDF, width = 14, height = 11)
print(p1); print(p2); print(p3)
invisible(grDevices::dev.off())

cat("\nwrote", OUT_CSV, "/", OUT_CLONE, "/", OUT_CELL, "/", OUT_PDF, "\nDIAG DONE\n")
