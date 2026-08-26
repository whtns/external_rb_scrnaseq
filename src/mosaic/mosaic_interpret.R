#!/usr/bin/env Rscript
# Issue #40 (interpretation add-on). The UNBIASED companion to mosaic_supplement.R.
# Where the supplement projects communities onto 4 hand-picked references
# (hypoxia + cell-cycle) and asks "which community matches these?", this script
# lets the data speak first:
#   P1  community x HALLMARK-50 ssGSEA NES         -> data-driven identity of every
#                                                     community, incl. the large
#                                                     pan-tumor #1 our 4 refs miss
#   P2  top marker genes per community (from the   -> read each community straight
#       representative program gene spectra)          from its own genes
#   P3  community x Seurat-cluster mean usage       -> does a program == a cluster?
#       (per sample; clusters are sample-specific)
#   P4  usage-on-UMAP for the pan-tumor communities -> where each recurrent program
#       (>=5/7 tumors), sample x community grid         is active in each tumor
#
# Deliverables: results/fig_mosaic_nmf_interpret.pdf,
#               results/table_mosaic_nmf_community_genes.csv,
#               docs/mosaic_nmf_interpret.md
#
# Usage: Rscript src/mosaic/mosaic_interpret.R <hallmark_ssgsea_dir> <usage_dir> <integrate_dir> <SRX> [<SRX> ...]
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(Seurat); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: mosaic_interpret.R <hallmark_ssgsea_dir> <usage_dir> <integrate_dir> <SRX> [...]")
ssgsea_dir <- args[[1]]; usage_dir <- args[[2]]; integ_dir <- args[[3]]; samples <- args[-(1:3)]

## ---- shared: parse a mosaicmpi ssGSEA representative_program_nes.txt ----------
# Format: header rows keyed community/dataset/k/program/Term/"" then one row per
# reference set (key = set name), columns = representative programs. We aggregate
# programs -> community by median NES. (Same shape mosaic_supplement.R parses.)
parse_nes_by_community <- function(nes_path) {
  raw <- strsplit(readLines(nes_path), "\t", fixed = TRUE)
  key <- vapply(raw, `[`, character(1), 1)
  hdr <- c("community", "dataset", "k", "program", "Term", "")
  comm_of_col <- raw[[which(key == "community")[1]]][-1]
  data_rows   <- which(!key %in% hdr & nzchar(key))
  lapply(data_rows, function(i) {
    v <- suppressWarnings(as.numeric(raw[[i]][-1]))
    tibble(reference = key[i], community = comm_of_col, nes = v)
  }) |> bind_rows() |> filter(!is.na(nes)) |>
    group_by(reference, community) |>
    summarise(nes = median(nes), .groups = "drop")
}

comm_num <- function(x) factor(x, levels = as.character(sort(unique(as.integer(x)))))

## ---- recurrence: which communities are pan-tumor (>=5/7 samples) --------------
pc <- fread(file.path(integ_dir, "program_communities.txt"), header = FALSE, sep = "\t")
pc_ds <- sub("\\|.*$", "", pc[[1]]); pc_comm <- as.character(pc[[2]])
recur <- tibble(community = pc_comm, dataset = pc_ds) |>
  distinct() |> count(community, name = "n_samples")
pan <- recur |> filter(n_samples >= 5) |> pull(community)
pan <- as.character(sort(as.integer(pan)))
message("pan-tumor communities (>=5 samples): ", paste(pan, collapse = ", "))

## ---- P1: community x HALLMARK NES --------------------------------------------
hm <- parse_nes_by_community(file.path(ssgsea_dir, "representative_program_nes.txt"))
# keep the Hallmark sets with real signal somewhere (|NES| >= 0.15 in >=1 community)
keep_sets <- hm |> group_by(reference) |>
  summarise(mx = max(abs(nes), na.rm = TRUE), .groups = "drop") |>
  filter(mx >= 0.15) |> pull(reference)
hm_f <- hm |> filter(reference %in% keep_sets) |>
  mutate(community = comm_num(community),
         reference = str_remove(reference, "^HALLMARK_"))
# order references by the community they peak in, for a readable diagonal-ish block
ord <- hm_f |> group_by(reference) |> slice_max(nes, n = 1, with_ties = FALSE) |>
  arrange(as.integer(as.character(community)), desc(nes)) |> pull(reference)
hm_f <- mutate(hm_f, reference = factor(reference, levels = rev(unique(ord))))
p1 <- ggplot(hm_f, aes(community, reference, fill = nes)) +
  geom_tile(color = "grey92") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "P1. Community identity vs MSigDB HALLMARK (ssGSEA NES, data-driven)",
       subtitle = "all 50 Hallmark sets scored; sets with |NES|>=0.15 shown",
       x = "mosaicMPI community", y = NULL, fill = "NES") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank())

## ---- P2: top genes per community (representative program spectra) ------------
rp <- fread(file.path(integ_dir, "representative_programs.txt"), header = FALSE, sep = "\t")
rp_comm <- as.character(unlist(rp[1, -1]))                 # community per column (row 1)
# header rows are Community/dataset/k/Program; genes are every row after them
meta_keys <- c("Community", "dataset", "k", "Program")
gene_rows <- which(!(rp[[1]] %in% meta_keys))
genes <- rp[[1]][gene_rows]
mat <- as.matrix(rp[gene_rows, -1]); storage.mode(mat) <- "numeric"
rownames(mat) <- genes
# z-score each program (column) so high-magnitude spectra don't dominate, then
# average per community -> relative gene importance.
matz <- scale(mat)
matz[!is.finite(matz)] <- 0
comm_levels <- as.character(sort(unique(as.integer(rp_comm))))
commz <- sapply(comm_levels, function(cc) rowMeans(matz[, rp_comm == cc, drop = FALSE], na.rm = TRUE))
rownames(commz) <- genes
# top-20 genes per community -> table
topN <- 20L
gene_tbl <- lapply(comm_levels, function(cc) {
  v <- sort(commz[, cc], decreasing = TRUE)
  tibble(community = cc, rank = seq_len(topN), gene = names(v)[seq_len(topN)],
         mean_z = round(as.numeric(v[seq_len(topN)]), 3))
}) |> bind_rows()
write_csv(gene_tbl, "results/table_mosaic_nmf_community_genes.csv")
cat("wrote results/table_mosaic_nmf_community_genes.csv\n")
# heatmap: top-5 genes per community (union) x community
top5 <- gene_tbl |> filter(rank <= 5) |> pull(gene) |> unique()
p2df <- as.data.frame(commz[top5, , drop = FALSE]) |>
  tibble::rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "community", values_to = "z") |>
  mutate(community = comm_num(community), gene = factor(gene, levels = rev(top5)))
p2 <- ggplot(p2df, aes(community, gene, fill = z)) +
  geom_tile(color = "grey92") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "P2. Top marker genes per community (representative program spectra, column-z)",
       x = "mosaicMPI community", y = NULL, fill = "mean z") +
  theme_minimal(base_size = 7) + theme(panel.grid = element_blank())

## ---- P3/P4: per-cell join to Seurat clusters + UMAP --------------------------
pick <- function(cands, have) { h <- intersect(cands, have); if (length(h)) h[1] else NA_character_ }
xcluster <- list(); umap_rows <- list()
for (srx in samples) {
  uf <- file.path(usage_dir, paste0("community_usage_", srx, ".csv"))
  sf <- file.path("output/seurat", paste0(srx, "_filtered_seu.rds"))
  if (!file.exists(uf) || !file.exists(sf)) { message("!! missing inputs for ", srx); next }
  usage <- read_csv(uf, show_col_types = FALSE)
  comm_cols <- setdiff(colnames(usage), "barcode")
  seu <- readRDS(sf)
  md <- seu@meta.data; md$barcode <- rownames(md)
  clcol <- pick(c("clusters", "seurat_clusters", "cell_type",
                  "SCT_snn_res.0.6", "SCT_snn_res.0.8"), colnames(md))
  if (is.na(clcol)) {                                   # any *_snn_res.* column
    snn <- grep("_snn_res", colnames(md), value = TRUE)
    clcol <- if (length(snn)) snn[1] else NA_character_
  }
  if (is.na(clcol)) { md$.cl <- as.character(Seurat::Idents(seu)); clcol <- ".cl" }
  # UMAP embedding (prefer a reduction named like umap)
  reds <- Seurat::Reductions(seu)
  ured <- reds[str_detect(tolower(reds), "umap")]; ured <- if (length(ured)) ured[1] else reds[1]
  emb <- as.data.frame(Seurat::Embeddings(seu, ured))[, 1:2]; colnames(emb) <- c("UMAP1", "UMAP2")
  emb$barcode <- rownames(emb)
  j <- usage |> inner_join(md[, c("barcode", clcol)], by = "barcode") |>
    inner_join(emb, by = "barcode")
  # P3: usage per (community, cluster). log1p tames mosaicMPI's un-normalized
  # usage -- a single degenerate cell can hit ~1e12, which destroys a raw mean and
  # pins any shared linear color scale so every other tile goes black. Mean of
  # log1p is robust to that tail; per-community relative scaling happens below.
  x <- j |> pivot_longer(all_of(comm_cols), names_to = "community", values_to = "usage") |>
    group_by(community, cluster = .data[[clcol]]) |>
    summarise(logu = mean(log1p(pmax(usage, 0)), na.rm = TRUE), .groups = "drop") |>
    mutate(sample = srx)
  xcluster[[srx]] <- x
  # P4: usage-on-UMAP for pan-tumor communities only (downsample cells so the
  # vector PDF stays small/fast -- the density pattern is unchanged).
  jj <- if (nrow(j) > 5000) j[sample.int(nrow(j), 5000), ] else j
  u <- jj |> select(barcode, UMAP1, UMAP2, all_of(intersect(pan, comm_cols))) |>
    pivot_longer(-c(barcode, UMAP1, UMAP2), names_to = "community", values_to = "usage") |>
    group_by(community) |> mutate(usage = usage / max(usage[is.finite(usage)], na.rm = TRUE)) |>
    ungroup() |> mutate(sample = srx)
  umap_rows[[srx]] <- u
  rm(seu, j); gc()
}
xcl <- bind_rows(xcluster); um <- bind_rows(umap_rows)
# Drop communities absent in a sample (all-NA -> non-finite mean, would render as
# solid gray columns), then min-max scale each community WITHIN each tumor to
# [0,1] so its cluster pattern is visible regardless of the large usage-scale
# differences between communities.
xcl <- xcl |> filter(is.finite(logu)) |>
  group_by(sample, community) |>
  mutate(rel = { r <- range(logu); if (diff(r) > 0) (logu - r[1]) / diff(r) else 0 }) |>
  ungroup()

p3 <- if (nrow(xcl)) ggplot(xcl, aes(comm_num(community), factor(cluster), fill = rel)) +
  geom_tile(color = "grey92") + facet_wrap(~ sample, scales = "free_y", ncol = 2) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey92") +
  labs(title = "P3. Community usage per Seurat cluster (per-community relative, within tumor)",
       subtitle = "mean log1p(usage) per cluster, min-max scaled within each community; communities absent in a tumor are omitted",
       x = "mosaicMPI community", y = "Seurat cluster", fill = "relative\nusage") +
  theme_minimal(base_size = 7) + theme(panel.grid = element_blank()) else NULL

p4 <- if (nrow(um)) ggplot(um, aes(UMAP1, UMAP2, color = usage)) +
  geom_point(size = 0.15, stroke = 0) +
  facet_grid(sample ~ paste0("c", community), switch = "y") +
  scale_color_viridis_c(option = "viridis") +
  labs(title = "P4. Pan-tumor community usage on UMAP (per-community max-scaled)",
       color = "usage") +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(size = 7),
        legend.position = "right",
        plot.title = element_text(size = 10, face = "bold")) else NULL

## ---- write multi-page PDF ----------------------------------------------------
dir.create("results", showWarnings = FALSE)
pdf("results/fig_mosaic_nmf_interpret.pdf", width = 11, height = 9)
print(p1); print(p2)
if (!is.null(p3)) print(p3)
if (!is.null(p4)) print(p4)
invisible(dev.off())
cat("wrote results/fig_mosaic_nmf_interpret.pdf\n")

## ---- methods md --------------------------------------------------------------
top_hallmark <- hm |> group_by(community) |> slice_max(nes, n = 1, with_ties = FALSE) |>
  transmute(community, top_hallmark = str_remove(reference, "^HALLMARK_"), nes = round(nes, 3)) |>
  arrange(as.integer(community))
dir.create("docs", showWarnings = FALSE)
writeLines(c(
  "# mosaicMPI NMF -- unbiased interpretation (issue #40 add-on)",
  "",
  "## Purpose",
  "Companion to the supplement: instead of projecting communities onto our 4 chosen",
  "hypoxia/cell-cycle references, characterise every community data-driven -- against",
  "the full 50-set MSigDB HALLMARK library (P1), from its own top marker genes (P2),",
  "against our Seurat clusters (P3), and spatially on each tumor's UMAP (P4).",
  "",
  sprintf("Pan-tumor communities (present in >=5/7 tumors): %s.", paste(pan, collapse = ", ")),
  "",
  "## Top HALLMARK match per community",
  paste0("| community | top HALLMARK | NES |"), "|---|---|---|",
  apply(top_hallmark, 1, function(r) sprintf("| %s | %s | %s |", r[["community"]], r[["top_hallmark"]], r[["nes"]])),
  "",
  "See results/fig_mosaic_nmf_interpret.pdf and results/table_mosaic_nmf_community_genes.csv."
), "docs/mosaic_nmf_interpret.md")
cat("wrote docs/mosaic_nmf_interpret.md\n")
cat("MOSAIC INTERPRET DONE\n")
