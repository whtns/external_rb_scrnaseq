# Diagnose "many very small clusters" in the low-hypoxia collages.
#
# Hypothesis: the collages display SCT_snn_res.* columns, but the low object is
# re-clustered by subset_seu_by_expression() on the GENE assay only
# (gene_snn_res.*). The SCT_snn_res.* columns are therefore INHERITED from the
# full pre-split object and merely subset to the low-hypoxia cells -- so any
# cluster whose cells were mostly high-hypoxia collapses to a tiny remnant.
#
# Test: for each paper_retained sample, report cluster-size stats for every
# *_snn_res.* column in the low object, separating SCT_ (inherited) from gene_
# (re-clustered post-split). If SCT_ is fragmented (many <20-cell clusters) while
# gene_ is clean, the small clusters are an artifact of not re-clustering SCT
# after the split -- NOT of the split's own (gene-assay) clustering.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr); library(readr); library(purrr)
})

paper <- c("SRX10264523","SRX10264526","SRX11133592","SRX11133593","SRX11133594")

size_stats <- function(meta, col, sample_id, object) {
  v <- as.character(meta[[col]])
  v <- v[!is.na(v)]
  tab <- sort(table(v), decreasing = TRUE)
  tibble(object = object, sample_id = sample_id, column = col,
         assay = str_extract(col, "^[^_]+"),
         res = as.numeric(str_extract(col, "[0-9.]+$")),
         n_cells = length(v), n_clusters = length(tab),
         min = min(tab), median = as.numeric(median(tab)),
         n_lt20 = sum(tab < 20), n_lt10 = sum(tab < 10), n_lt5 = sum(tab < 5))
}

all_rows <- list()
for (sid in paper) {
  low <- sprintf("output/seurat/%s_hypoxia_low_seu.rds", sid)
  if (!file.exists(low)) { cat("MISSING low:", low, "\n"); next }
  s <- readRDS(low)
  meta <- s@meta.data
  cols <- grep("_snn_res\\.[0-9.]+$", colnames(meta), value = TRUE)
  cat("\n==", sid, "-- low object:", ncol(s), "cells; snn cols:",
      paste(cols, collapse = ", "), "\n")
  for (col in cols) all_rows[[paste(sid, col)]] <- size_stats(meta, col, sid, "low")
  rm(s, meta); gc()
}

res_tbl <- bind_rows(all_rows) %>% arrange(sample_id, assay, res)
out <- "results/hypoxia_cluster_split/low_hypoxia_cluster_size_diag.csv"
write_csv(res_tbl, out)

cat("\n\n========= SCT (inherited, shown in collages) vs gene (re-clustered) =========\n")
res_tbl %>%
  filter(assay %in% c("SCT","gene")) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("\nWrote", out, "\n")
