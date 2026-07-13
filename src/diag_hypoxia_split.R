#!/usr/bin/env Rscript
# Validate that the cluster-based high partition is genuinely hypoxic:
# compare high vs low on hypoxia_score and on marker-gene mean expression.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
labeled <- if (length(args) >= 1) args[[1]] else
  "output/seurat/SRX10031194_seu_hypoxia_labeled.rds"

seu <- readRDS(labeled)
DefaultAssay(seu) <- "gene"
part <- seu$hypoxia_partition

cat("== hypoxia_score by partition ==\n")
print(tapply(seu$hypoxia_score, part, function(x) round(mean(x, na.rm = TRUE), 3)))
cat("\n== n cells ==\n"); print(table(part))

genes <- c("ZFAS1","GAS5","SNHG5","SNHG6","SNHG7","SNHG8","SNHG15",
           "VEGFA","NDRG1","BNIP3","SLC2A1","CA9","PGK1","LDHA","ENO1","P4HA1")
genes <- genes[genes %in% rownames(seu)]
expr <- GetAssayData(seu, assay = "gene", layer = "data")[genes, , drop = FALSE]

cat("\n== mean log-norm expression: high vs low (sorted by high) ==\n")
hi <- Matrix::rowMeans(expr[, part == "high", drop = FALSE])
lo <- Matrix::rowMeans(expr[, part == "low",  drop = FALSE])
df <- data.frame(gene = genes, high = round(hi, 3), low = round(lo, 3),
                 diff = round(hi - lo, 3))
df <- df[order(-df$diff), ]
print(df, row.names = FALSE)

cat("\n== module score of marker set, high vs low ==\n")
seu <- AddModuleScore(seu, features = list(hyp = genes), name = "hypset", nbin = 12, ctrl = 50)
print(tapply(seu$hypset1, part, function(x) round(mean(x), 3)))
cat("\nDONE\n")
