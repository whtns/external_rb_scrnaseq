suppressMessages({library(targets); library(dplyr); library(stringr)}); tar_config_set(store="_targets_r431")
lcc <- tar_read_raw("large_clone_comparisons")
md  <- readr::read_tsv("data/metadata.tsv", show_col_types=FALSE) |> select(Run, Experiment) |> distinct()
srx2srr <- setNames(md$Run, md$Experiment)
samples <- c("SRX11133592","SRX11133593","SRX11133594","SRX10264523","SRX10264526","SRX10831287")
for (s in samples) {
  srr <- srx2srr[[s]]
  keys <- names(lcc)[str_detect(names(lcc), fixed(srr %||% "NOPE"))]
  cmp  <- unlist(lapply(lcc[keys], names))
  cat(sprintf("%-12s -> SRR=%-12s keys=%-45s\n", s, srr %||% "NONE", paste(keys, collapse=",")))
  cat("             comparisons:", paste(cmp, collapse=" | "), "\n")
  cat("             1q+ hits:", sum(str_detect(cmp, fixed("1q+"))), "| 16q- hits:", sum(str_detect(cmp, fixed("16q-"))), "\n")
}
