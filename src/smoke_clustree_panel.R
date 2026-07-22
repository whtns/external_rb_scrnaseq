# Confirm the clustree collage panel (#38) renders from a multi-resolution sweep.
# The low-hypoxia objects carry persisted SCT_snn_res.0.2..1.0, so
# .build_clustree_panel() should return a ggplot; save it to prove it renders.
suppressPackageStartupMessages({
  source("packages.R")
  library(Seurat); library(ggplot2)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})

p   <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
seu <- readRDS(p)
cols <- grep("^SCT_snn_res\\.", colnames(seu@meta.data), value = TRUE)
cat("resolution columns present:", paste(cols, collapse = ", "), "\n")

pnl <- numbatHelpers:::.build_clustree_panel(seu, assay = "SCT")
cat("clustree panel class:", paste(class(pnl), collapse = "/"), "\n")
stopifnot(inherits(pnl, "ggplot"))

out <- "tmp/smoke_clustree_panel.pdf"
fs::dir_create("tmp")
ggsave(out, pnl, width = 6, height = 6)
sz <- file.info(out)$size
cat(sprintf("wrote %s | %s bytes | %s\n", out, format(sz, big.mark = ","),
            if (sz > 10000) "OK" else "!! TOO SMALL"))
cat("CLUSTREE SMOKE DONE\n")
