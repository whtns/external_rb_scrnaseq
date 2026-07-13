#!/usr/bin/env Rscript
# Audit the abbreviation annotations in data/cluster_dictionary.tsv against the
# actual gene_snn_res.0.2 marker genes of each unfiltered Seurat object.
#
# Primary question: do all MALAT1-annotated clusters have MALAT1 as a (top)
# marker gene? More generally, produce, for every (sample, cluster), the top
# marker genes so we can both debug and later automate the dictionary.
#
# Output: doc/cluster_marker_audit.tsv  (one row per sample/cluster)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

devtools::load_all("/project2/cobrinik_1090/rpkgs/seuratTools", quiet = TRUE)

dict <- read_tsv("data/cluster_dictionary.tsv", show_col_types = FALSE)
samples <- sort(unique(dict$sample_id))

# how many top markers to record per cluster
TOP_N <- 10
RES_COL <- "gene_snn_res.0.2"

get_markers <- function(seu) {
  # Prefer stashed presto markers; otherwise compute them.
  m <- seu@misc$markers[[RES_COL]]$presto
  if (is.null(m)) {
    if (!RES_COL %in% colnames(seu@meta.data)) return(NULL)
    DefaultAssay(seu) <- "gene"
    if (inherits(seu[["gene"]], "Assay5")) seu[["gene"]] <- JoinLayers(seu[["gene"]])
    m <- presto::wilcoxauc(seu, RES_COL, seurat_assay = "gene") %>%
      dplyr::group_by(group) %>%
      dplyr::filter(padj < 0.5) %>%
      dplyr::arrange(group, desc(logFC)) %>%
      dplyr::select(Gene.Name = feature,
                    Average.Log.Fold.Change = logFC,
                    Adjusted.pvalue = padj,
                    avgExpr,
                    Cluster = group) %>%
      dplyr::ungroup()
  }
  m
}

rows <- list()
for (s in samples) {
  seu_path <- sprintf("output/seurat/%s_seu.rds", s)
  if (!file.exists(seu_path)) { message("skip (no seu): ", s); next }
  message("=== ", s, " ===")
  seu <- readRDS(seu_path)
  m <- get_markers(seu)
  if (is.null(m)) { message("  no ", RES_COL, " / markers for ", s); rm(seu); gc(); next }

  m <- m %>%
    mutate(Cluster = as.character(Cluster)) %>%
    arrange(Cluster, desc(Average.Log.Fold.Change))

  sdict <- dict %>% filter(sample_id == s) %>%
    mutate(gene_snn_res.0.2 = as.character(gene_snn_res.0.2))

  clusters <- sort(unique(m$Cluster))
  for (cl in clusters) {
    cm <- m %>% filter(Cluster == cl)
    top <- head(cm, TOP_N)
    # rank of MALAT1 within this cluster's marker list (by logFC)
    malat_idx <- which(cm$Gene.Name == "MALAT1")
    malat_rank <- if (length(malat_idx)) malat_idx[1] else NA_integer_
    malat_lfc  <- if (length(malat_idx)) cm$Average.Log.Fold.Change[malat_idx[1]] else NA_real_
    malat_padj <- if (length(malat_idx)) cm$Adjusted.pvalue[malat_idx[1]] else NA_real_

    abbr_row <- sdict %>% filter(gene_snn_res.0.2 == cl)
    abbr <- if (nrow(abbr_row)) abbr_row$abbreviation[1] else NA_character_
    rmv  <- if (nrow(abbr_row)) abbr_row$remove[1] else NA

    rows[[length(rows) + 1]] <- tibble(
      sample_id        = s,
      cluster          = cl,
      current_abbrev   = abbr,
      current_remove   = rmv,
      top1_gene        = top$Gene.Name[1],
      top1_lfc         = top$Average.Log.Fold.Change[1],
      top_genes        = paste(top$Gene.Name, collapse = ", "),
      malat1_rank      = malat_rank,
      malat1_lfc       = malat_lfc,
      malat1_padj      = malat_padj
    )
  }
  rm(seu, m); gc()
}

audit <- bind_rows(rows)
out <- "doc/cluster_marker_audit.tsv"
write_tsv(audit, out)
message("wrote ", out, " (", nrow(audit), " rows)")

# ---- focused MALAT1 report ----
message("\n================ MALAT1-annotated clusters ================")
malat_clusters <- audit %>% filter(current_abbrev == "MALAT1")
print(as.data.frame(malat_clusters %>%
  select(sample_id, cluster, current_remove, top1_gene, top1_lfc,
         malat1_rank, malat1_lfc)))

message("\nDo all MALAT1-annotated clusters have MALAT1 as top marker? ",
        all(malat_clusters$top1_gene == "MALAT1", na.rm = TRUE))
message("Do all MALAT1-annotated clusters have MALAT1 present as any marker? ",
        all(!is.na(malat_clusters$malat1_rank)))
