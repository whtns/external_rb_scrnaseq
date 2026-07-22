# Issue #31: how much of each 6p call actually rests on HLA / MHC genes?
#
# Before rerunning anything with a masked gtf, measure the exposure from the
# numbat output we already have. Masking can only change a call in proportion to
# what the masked features contribute, and numbat's own output decomposes each
# segment's evidence into LLR_y (expression) and LLR_x (allele) -- which matters
# here because gtf masking touches ONLY the expression arm:
#
#   run_numbat() line 29-31:  count_mat = count_mat[intersect(gtf$gene, ...), ]
#     -> dropping HLA rows from the gtf DOES remove them from the expression model
#   run_numbat() line 26:     df_allele = annotate_genes(df_allele, gtf)
#     -> a left_join; SNPs that hit no gtf gene get gene = NA and are KEPT, and
#        get_bulk() then explicitly keeps them: filter(lambda_ref != 0 | is.na(gene))
#     -> MHC SNPs still drive the allele arm after gtf masking
#
# So the report below splits the two arms deliberately. A segment whose evidence
# is mostly LLR_x cannot be rescued by a gtf-only mask no matter how many HLA
# genes it contains.
#
# Writes results/hla_6p_contribution.csv (one row per non-neutral chr6p segment).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr); library(tidyr)
  library(numbat)
})

NUMBAT_DIR <- "output/numbat_sridhar"

# HLA-prefixed genes, and the wider MHC window. The window is the thing the issue
# asks about in point 3: MICA/MICB/TAP1/TAP2/complement share the polymorphism and
# interferon-response problems without carrying an "HLA-" prefix.
gtf <- as.data.frame(numbat::gtf_hg38)
hla_genes <- gtf$gene[grepl("^HLA-", gtf$gene) & gtf$CHROM == 6]
MHC_START <- 28.5e6
MHC_END   <- 33.5e6
mhc_genes <- gtf$gene[gtf$CHROM == 6 & gtf$gene_end > MHC_START & gtf$gene_start < MHC_END]
cat("HLA- genes in gtf_hg38:", length(hla_genes), "\n")
cat("genes in MHC window", MHC_START/1e6, "-", MHC_END/1e6, "Mb:", length(mhc_genes), "\n")
cat("  non-HLA-prefixed among them:", length(setdiff(mhc_genes, hla_genes)), "\n")
cat("  e.g.:", paste(head(setdiff(mhc_genes, hla_genes), 15), collapse = ", "), "\n\n")

# chr6 p-arm ends at the centromere.
acen <- as.data.frame(numbat::acen_hg38)
cen6 <- acen[as.character(acen$CHROM) == "6", ]
P_END <- min(cen6$start)
cat("chr6 centromere start (p-arm end):", round(P_END / 1e6, 2), "Mb\n\n")

sample_dirs <- list.dirs(NUMBAT_DIR, recursive = FALSE)

one_sample <- function(d) {
  sid <- basename(d)
  bulk_f <- file.path(d, "bulk_clones_final.tsv.gz")
  if (!file.exists(bulk_f)) return(NULL)

  bulk <- tryCatch(
    readr::read_tsv(bulk_f, show_col_types = FALSE, progress = FALSE,
                    col_types = readr::cols(.default = readr::col_guess(),
                                            CHROM = readr::col_character(),
                                            seg = readr::col_character(),
                                            gene = readr::col_character())),
    error = function(e) { message("!! ", sid, ": ", conditionMessage(e)); NULL })
  if (is.null(bulk) || !all(c("CHROM", "seg", "cnv_state") %in% names(bulk))) return(NULL)

  # Non-neutral segments on the p-arm of chr6. `cnv_state_post` is the call numbat
  # reports; fall back to cnv_state if the column is absent.
  state_col <- if ("cnv_state_post" %in% names(bulk)) "cnv_state_post" else "cnv_state"

  b6 <- bulk %>%
    filter(CHROM == "6", !is.na(seg), .data[[state_col]] != "neu", seg_start < P_END)
  if (nrow(b6) == 0) return(NULL)

  b6 %>%
    group_by(seg) %>%
    summarise(
      sample_id  = sid,
      cnv_state  = first(.data[[state_col]]),
      seg_start  = first(seg_start),
      seg_end    = first(seg_end),
      LLR        = first(LLR),
      LLR_x      = first(LLR_x),      # allele arm -- gtf masking does NOT touch this
      LLR_y      = first(LLR_y),      # expression arm -- gtf masking DOES
      phi_mle    = first(phi_mle),
      # --- expression arm exposure -------------------------------------------
      n_genes_seg = n_distinct(gene[!is.na(gene)]),
      n_hla       = n_distinct(gene[gene %in% hla_genes]),
      n_mhc       = n_distinct(gene[gene %in% mhc_genes]),
      Y_total     = sum(Y_obs[!duplicated(gene) & !is.na(gene)], na.rm = TRUE),
      Y_hla       = sum(Y_obs[!duplicated(gene) & gene %in% hla_genes], na.rm = TRUE),
      Y_mhc       = sum(Y_obs[!duplicated(gene) & gene %in% mhc_genes], na.rm = TRUE),
      # --- allele arm exposure ------------------------------------------------
      n_snps_seg  = n_distinct(snp_id),
      n_snps_mhc  = n_distinct(snp_id[POS >= MHC_START & POS <= MHC_END]),
      DP_total    = sum(DP[!duplicated(snp_id)], na.rm = TRUE),
      DP_mhc      = sum(DP[!duplicated(snp_id) & POS >= MHC_START & POS <= MHC_END], na.rm = TRUE),
      .groups = "drop"
    )
}

res <- purrr::map(sample_dirs, one_sample) %>% purrr::compact() %>% bind_rows()

if (nrow(res) == 0) { cat("no non-neutral chr6p segments found\n"); quit(save = "no") }

res <- res %>%
  mutate(
    overlaps_mhc = seg_end > MHC_START & seg_start < MHC_END,
    frac_genes_hla = n_hla / pmax(n_genes_seg, 1),
    frac_genes_mhc = n_mhc / pmax(n_genes_seg, 1),
    frac_Y_hla     = Y_hla / pmax(Y_total, 1),
    frac_Y_mhc     = Y_mhc / pmax(Y_total, 1),
    frac_snps_mhc  = n_snps_mhc / pmax(n_snps_seg, 1),
    frac_DP_mhc    = DP_mhc / pmax(DP_total, 1),
    # Share of total evidence sitting in the arm a gtf mask can reach at all.
    frac_LLR_expr  = LLR_y / pmax(LLR_x + LLR_y, 1e-9)
  ) %>%
  arrange(desc(overlaps_mhc), desc(frac_Y_mhc))

readr::write_csv(res, "results/hla_6p_contribution.csv")
cat("wrote results/hla_6p_contribution.csv --", nrow(res), "segments across",
    n_distinct(res$sample_id), "samples\n\n")

cat("=== chr6p non-neutral segments that OVERLAP the MHC window ===\n")
ov <- res %>% filter(overlaps_mhc)
if (nrow(ov) == 0) {
  cat("NONE -- no called 6p segment overlaps 28.5-33.5 Mb, so HLA masking\n")
  cat("cannot change any current 6p call.\n")
} else {
  ov %>%
    transmute(sample_id, seg, cnv_state,
              Mb = sprintf("%.1f-%.1f", seg_start/1e6, seg_end/1e6),
              LLR = round(LLR, 1), LLR_x = round(LLR_x, 1), LLR_y = round(LLR_y, 1),
              pct_LLR_expr = round(100 * frac_LLR_expr),
              genes = n_genes_seg, hla = n_hla, mhc = n_mhc,
              pct_Y_hla = round(100 * frac_Y_hla, 1),
              pct_Y_mhc = round(100 * frac_Y_mhc, 1),
              pct_DP_mhc = round(100 * frac_DP_mhc, 1)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

cat("\n=== chr6p non-neutral segments that DO NOT reach the MHC ===\n")
nov <- res %>% filter(!overlaps_mhc)
cat(nrow(nov), "segments across", n_distinct(nov$sample_id),
    "samples -- unaffected by any HLA/MHC mask by construction\n")
if (nrow(nov) > 0) {
  nov %>%
    transmute(sample_id, seg, cnv_state,
              Mb = sprintf("%.1f-%.1f", seg_start/1e6, seg_end/1e6),
              LLR = round(LLR, 1)) %>%
    as.data.frame() %>% head(30) %>% print(row.names = FALSE)
}

cat("\n=== summary over MHC-overlapping segments ===\n")
if (nrow(ov) > 0) {
  cat("median % of segment expression counts from HLA- genes :",
      round(100 * median(ov$frac_Y_hla), 2), "%\n")
  cat("median % from the whole MHC window                    :",
      round(100 * median(ov$frac_Y_mhc), 2), "%\n")
  cat("median % of segment SNP depth inside the MHC window   :",
      round(100 * median(ov$frac_DP_mhc), 2), "%\n")
  cat("median % of segment LLR carried by the EXPRESSION arm :",
      round(100 * median(ov$frac_LLR_expr)), "%\n")
  cat("  (a gtf-only mask can move at most the expression arm)\n")
}

cat("\nDONE\n")
