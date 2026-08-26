#!/usr/bin/env Rscript
# Diagnose why plot_scna_two_clone_res_collages skips most FILTERED samples while the
# low-hypoxia twins succeed (issue #36 follow-up). Prints, per sample x scna, the
# config key resolution, the matched comparison, the retained clones, and the cell
# count in those clones -- i.e. every guard the function can exit on.
suppressPackageStartupMessages({
  source("packages.R")
  library(targets); library(Seurat); library(stringr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")

lcc <- tar_read(large_clone_comparisons)
retained <- tar_read(paper_retained_samples)

pairs <- list(
  c("SRX10264523", "1q"), c("SRX10264523", "2p"),
  c("SRX10264526", "1q"), c("SRX10264526", "2p"),
  c("SRX11133593", "1q"), c("SRX11133593", "16q"),
  c("SRX11133594", "1q"), c("SRX11133594", "16q")
)

for (pr in pairs) {
  srx <- pr[[1]]; scna <- pr[[2]]
  for (kind in c("filtered", "hypoxia_low")) {
    p <- if (kind == "filtered")
      sprintf("output/seurat/%s_filtered_seu.rds", srx)
    else
      sprintf("output/seurat/%s_hypoxia_low_seu.rds", srx)
    tag <- sprintf("%-12s %-4s %-11s", srx, scna, kind)
    if (!file.exists(p)) { cat(tag, "| NO FILE\n"); next }

    slug <- str_remove(fs::path_file(p), "_filtered_seu.*")
    key  <- if (slug %in% names(lcc)) slug else srx
    comps <- names(lcc[[key]])
    comp  <- comps[str_detect(comps, fixed(scna))]
    if (length(comp) == 0) {
      cat(tag, sprintf("| key=%s | NO %s COMPARISON (have: %s)\n",
                       key, scna, paste(comps, collapse=",")))
      next
    }
    rc <- comp |> str_extract("[0-9]+_v_[0-9]+") |> str_split("_v_", simplify = TRUE) |> as.vector()
    rc <- unique(rc[!is.na(rc) & rc != ""])

    seu <- readRDS(p)
    if (!"clone_opt" %in% colnames(seu@meta.data)) {
      cat(tag, sprintf("| key=%s comp=%s | NO clone_opt COLUMN\n", key, comp[1])); rm(seu); gc(); next
    }
    cl <- as.character(seu@meta.data$clone_opt)
    n_keep <- sum(cl %in% rc, na.rm = TRUE)
    cat(tag, sprintf("| key=%s comp=%s clones=%s | cells_in_clones=%d %s | clone_tab=%s\n",
                     key, comp[1], paste(rc, collapse="/"), n_keep,
                     if (n_keep < 20) "<< SKIP (<20)" else "ok",
                     paste(names(table(cl)), table(cl), sep=":", collapse=" ")))
    rm(seu); gc()
  }
}
cat("\nDIAG DONE\n")
