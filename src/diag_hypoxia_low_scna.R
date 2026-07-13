#!/usr/bin/env Rscript
# Diagnostic: why do the hypoxia_low cluster facet plots have no points, and
# why is `scna` absent/empty in metadata? Reports, per sample: cell count,
# whether `scna` exists, its value table, and NA counts for the facet plot's
# point aesthetics (S.Score, G2M.Score) plus clusters.
suppressPackageStartupMessages({
  library(Seurat)
})

files <- c(
  "output/seurat/SRX10831286_hypoxia_low_seu.rds",  # ABSENT from simplifications yaml
  "output/seurat/SRX11133589_hypoxia_low_seu.rds",  # PRESENT in yaml (1q+,16q-,17-,19p-)
  "output/seurat/SRX22868105_hypoxia_low_seu.rds",  # PRESENT in yaml (6p+,16q-)
  "output/seurat/SRX10831282_hypoxia_low_seu.rds"   # ABSENT from yaml
)

for (f in files) {
  cat("\n================================================================\n")
  cat(f, "\n")
  if (!file.exists(f)) { cat("  MISSING FILE\n"); next }
  seu <- readRDS(f)
  md  <- seu@meta.data
  cat("  n_cells:", ncol(seu), "\n")
  cat("  'scna' in meta.data:", "scna" %in% colnames(md), "\n")
  if ("scna" %in% colnames(md)) {
    cat("  scna value table (incl NA):\n")
    print(table(md$scna, useNA = "always"))
    cat("  scna: n_NA =", sum(is.na(md$scna)),
        " n_nonNA =", sum(!is.na(md$scna)),
        " n_blank('') =", sum(md$scna %in% ""), "\n")
  }
  for (col in c("S.Score", "G2M.Score", "clusters")) {
    if (col %in% colnames(md)) {
      v <- md[[col]]
      cat(sprintf("  %-10s: n_NA=%d  n_nonNA=%d  class=%s\n",
                  col, sum(is.na(v)), sum(!is.na(v)), class(v)[1]))
    } else {
      cat(sprintf("  %-10s: ABSENT\n", col))
    }
  }
  rm(seu, md); gc(FALSE)
}
cat("\nDONE\n")

# --- follow-up: confirm scores==0 in subset vs real scores in parent ---
cat("\n===== SCORE/SCNA: subset vs parent =====\n")
pairs <- list(
  c("output/seurat/SRX10831286_hypoxia_low_seu.rds", "output/seurat/SRX10831286_filtered_seu.rds"),
  c("output/seurat/SRX11133589_hypoxia_low_seu.rds", "output/seurat/SRX11133589_filtered_seu.rds")
)
for (pr in pairs) {
  for (f in pr) {
    if (!file.exists(f)) { cat(f, ": MISSING\n"); next }
    md <- readRDS(f)@meta.data
    ss <- if ("S.Score" %in% names(md)) md$S.Score else NA
    gg <- if ("G2M.Score" %in% names(md)) md$G2M.Score else NA
    sc <- if ("scna" %in% names(md)) md$scna else NA
    cat(sprintf("%s\n  S.Score:  %s\n  G2M.Score:%s\n  scna nonblank: %d/%d\n",
      f,
      paste(round(range(ss, na.rm=TRUE),4), collapse=" .. "),
      paste(round(range(gg, na.rm=TRUE),4), collapse=" .. "),
      sum(!is.na(sc) & sc != ""), length(sc)))
    rm(md); gc(FALSE)
  }
}
cat("\nDONE2\n")
