#!/usr/bin/env Rscript
# Issue #40, Stage G. Assemble the mosaicMPI NMF supplement from the fresh re-run:
#   1. results/fig_mosaic_nmf_supplement.pdf
#        A) gene-identity view  - representative-program x reference-set ssGSEA NES,
#           aggregated per community (which communities ARE hypoxia / cell-cycle).
#        B) per-cell view       - each community's per-cell usage correlated with
#           hypoxia_score / S.Score / G2M.Score, per sample (within-tumor match) and
#           pooled (cross-tumor recurrence).
#   2. results/table_mosaic_nmf_programs.csv  - per community summary.
#   3. docs/mosaic_nmf_supplement.md          - methods + headline results.
#
# Inputs (produced by earlier stages):
#   <ssgsea_dir>/representative_program_nes.txt   (mosaicmpi ssgsea -n ...)
#   <usage_dir>/community_usage_<SRX>.csv         (export_community_usage.py)
#   output/seurat/<SRX>_filtered_seu.rds          (per-cell hypoxia_score/S.Score/G2M.Score)
#
# Usage: Rscript src/mosaic/mosaic_supplement.R <ssgsea_dir> <usage_dir> <SRX> [<SRX> ...]
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(readr); library(stringr); library(Seurat)
  # numbatHelpers supplies .ensure_cc_scores() and add_hypoxia_score(), the SAME
  # scoring the pipeline uses -- so the per-cell scores we correlate against are
  # exactly "the scores we already use", computed on the mosaicMPI input cells.
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: mosaic_supplement.R <ssgsea_dir> <usage_dir> <SRX> [<SRX> ...]")
ssgsea_dir <- args[[1]]; usage_dir <- args[[2]]; samples <- args[-(1:2)]

score_cols <- c(hypoxia = "hypoxia_score", S = "S.Score", G2M = "G2M.Score")

# Add S.Score/G2M.Score/hypoxia_score to a filtered_seu that lacks them, so the
# per-cell correlation panel is populated (the filtered_seu objects frequently
# carry none of these -- only the hypoxia-split objects do).
ensure_scores <- function(seu) {
  seu <- numbatHelpers:::.ensure_cc_scores(seu)
  if (!"hypoxia_score" %in% colnames(seu@meta.data)) {
    seu <- tryCatch(add_hypoxia_score(seu),
                    error = function(e) { message("add_hypoxia_score failed: ",
                                                   conditionMessage(e)); seu })
  }
  seu
}

## ---- A. representative-program ssGSEA NES -> reference x community ------------
nes_path <- file.path(ssgsea_dir, "representative_program_nes.txt")
if (!file.exists(nes_path)) stop("missing ssGSEA output: ", nes_path)
raw <- strsplit(readLines(nes_path), "\t", fixed = TRUE)
key <- vapply(raw, `[`, character(1), 1)
hdr <- c("community", "dataset", "k", "program", "Term", "")
comm_of_col <- raw[[which(key == "community")[1]]][-1]           # community label per column
data_rows   <- which(!key %in% hdr & nzchar(key))
ref_x_comm <- lapply(data_rows, function(i) {
  v <- suppressWarnings(as.numeric(raw[[i]][-1]))
  tibble(reference = key[i], community = comm_of_col, nes = v)
}) %>% bind_rows() %>%
  filter(!is.na(nes)) %>%
  group_by(reference, community) %>%
  summarise(nes = median(nes), .groups = "drop")            # aggregate programs -> community

# Headline references first, then the rest.
headline <- c("HALLMARK_HYPOXIA", "TIROSH_S_PHASE", "TIROSH_G2M", "GIOTTI_CELL_CYCLE")
ref_levels <- c(intersect(headline, ref_x_comm$reference),
                setdiff(sort(unique(ref_x_comm$reference)), headline))
ref_x_comm <- mutate(ref_x_comm,
                     reference = factor(reference, levels = rev(ref_levels)),
                     community = factor(community, levels = sort(unique(as.integer(community)))))

panelA <- ggplot(filter(ref_x_comm, as.character(reference) %in% headline),
                 aes(community, reference, fill = nes)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "A. Representative program ssGSEA NES vs reference programs",
       x = "mosaicMPI community", y = NULL, fill = "NES") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank())

## ---- B. per-cell community usage vs the scores we already use ----------------
cor_rows <- list()
for (srx in samples) {
  uf <- file.path(usage_dir, paste0("community_usage_", srx, ".csv"))
  sf <- file.path("output/seurat", paste0(srx, "_filtered_seu.rds"))
  if (!file.exists(uf)) { message("!! no usage for ", srx, "; skipping"); next }
  usage <- readr::read_csv(uf, show_col_types = FALSE)
  comm_cols <- setdiff(colnames(usage), "barcode")
  seu <- readRDS(sf)
  seu <- ensure_scores(seu)
  md <- seu@meta.data
  md$barcode <- rownames(md)
  present <- score_cols[score_cols %in% colnames(md)]  # keep names (axis labels)
  # Hard barcode-identity assertion: every usage row must join to a Seurat cell.
  n_join <- length(intersect(usage$barcode, md$barcode))
  if (n_join < nrow(usage))
    stop(sprintf("barcode mismatch for %s: %d/%d usage rows join to Seurat cells",
                 srx, n_join, nrow(usage)))
  j <- dplyr::inner_join(usage, md[, c("barcode", present)], by = "barcode")
  for (cc in comm_cols) for (ax in names(present)) {
    v <- j[[cc]]; s <- j[[present[[ax]]]]
    ok <- is.finite(v) & is.finite(s)
    r <- if (sum(ok) > 10 && sd(v[ok]) > 0 && sd(s[ok]) > 0)
      suppressWarnings(cor(v[ok], s[ok], method = "spearman")) else NA_real_
    cor_rows[[length(cor_rows) + 1]] <-
      tibble(sample = srx, community = cc, axis = ax, r = r, n = sum(ok))
  }
  rm(seu); gc()
}
cor_tbl <- bind_rows(cor_rows)

panelB <- if (nrow(cor_tbl)) {
  ggplot(cor_tbl, aes(factor(community), sample, fill = r)) +
    geom_tile(color = "grey90") +
    facet_wrap(~ axis, ncol = 1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "B. Per-cell community usage vs pipeline scores (Spearman)",
         x = "mosaicMPI community", y = NULL, fill = "rho") +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
} else patchwork::plot_spacer()

fig <- (panelA / panelB) + plot_layout(heights = c(1, 2)) +
  plot_annotation(title = "mosaicMPI NMF programs vs hypoxia / Tirosh cell-cycle references")
dir.create("results", showWarnings = FALSE)
ggsave("results/fig_mosaic_nmf_supplement.pdf", fig, width = 12, height = 12)
cat("wrote results/fig_mosaic_nmf_supplement.pdf\n")

## ---- table: per community summary --------------------------------------------
best_ref <- ref_x_comm %>% group_by(community) %>%
  slice_max(nes, n = 1, with_ties = FALSE) %>%
  transmute(community, best_reference = as.character(reference), best_reference_nes = nes)
best_cor <- if (nrow(cor_tbl)) cor_tbl %>%
  group_by(community) %>% slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  transmute(community, top_axis = axis, top_axis_rho = r, top_axis_sample = sample) else
  tibble(community = character())
tbl <- best_ref %>%
  mutate(community = as.character(community)) %>%
  full_join(mutate(best_cor, community = as.character(community)), by = "community") %>%
  arrange(suppressWarnings(as.integer(community)))
readr::write_csv(tbl, "results/table_mosaic_nmf_programs.csv")
cat("wrote results/table_mosaic_nmf_programs.csv (", nrow(tbl), " communities)\n", sep = "")

## ---- methods md --------------------------------------------------------------
hy_comm <- best_ref %>% filter(best_reference == "HALLMARK_HYPOXIA") %>% pull(community)
g2m_comm <- best_ref %>% filter(best_reference == "TIROSH_G2M") %>% pull(community)
dir.create("docs", showWarnings = FALSE)
writeLines(c(
  "# mosaicMPI NMF gene-program supplement (issue #40)",
  "",
  "## Methods",
  sprintf("Consensus NMF (cNMF, via mosaicMPI) was run on the %d paper-retained tumors (%s).",
          length(samples), paste(samples, collapse = ", ")),
  "Raw gene counts per tumor were factorized over a rank sweep, programs were consensus-",
  "aggregated (postprocess), and programs were integrated across tumors into a community",
  "network (single, SCNA-agnostic integration). Recovered programs were compared to the",
  "signatures the main analysis uses: HALLMARK_HYPOXIA, Tirosh S-phase and G2/M markers",
  "(Seurat::cc.genes.updated.2019), and the giotti cell-cycle gene sets. Comparison used",
  "(i) ssGSEA NES of each representative program against those reference sets, aggregated",
  "per community, and (ii) Spearman correlation of each community's per-cell usage",
  "(normalize=False) with the per-cell hypoxia_score / S.Score / G2M.Score.",
  "",
  "## Result",
  sprintf("Community/communities best matching HALLMARK_HYPOXIA: %s.",
          if (length(hy_comm)) paste(hy_comm, collapse = ", ") else "none"),
  sprintf("Community/communities best matching TIROSH_G2M: %s.",
          if (length(g2m_comm)) paste(g2m_comm, collapse = ", ") else "none"),
  "See results/fig_mosaic_nmf_supplement.pdf and results/table_mosaic_nmf_programs.csv."
), "docs/mosaic_nmf_supplement.md")
cat("wrote docs/mosaic_nmf_supplement.md\n")
cat("MOSAIC SUPPLEMENT DONE\n")
