# Smoke-test find_diffex_bw_clones_for_each_cluster() on ONE low-hypoxia object
# before the full diffex rebuild. Confirms the two guarded failure modes are gone:
#   (1) cis "`x` must be a vector" (missing `clusters` col -> derived from group.by)
#   (2) the NA-path crash in the downstream tabulator (.drop_missing_paths)
# and that a real per-cluster/clone diffex CSV is produced on the CURRENT
# (post-07-17 fresh-PCA) low-object clusters.
suppressPackageStartupMessages({
  source("packages.R"); library(Seurat); library(targets)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")

lcc <- tar_read(large_clone_comparisons)
nbf <- tar_read(numbat_rds_files)
co  <- tar_read(cluster_orders)
low <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
stopifnot(file.exists(low))
cat("comparisons for SRX10264523:",
    paste(names(lcc[["SRX10264523"]]), collapse = ", "), "\n")

paths <- c()
for (loc in c("all", "cis", "out_of_segment")) {
  cat("\n=== location:", loc, "===\n")
  out <- tryCatch(
    find_diffex_bw_clones_for_each_cluster(low, nbf, lcc, co, location = loc),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("returned:", out, "\n")
  if (!is.na(out) && file.exists(out)) {
    d <- readr::read_csv(out, show_col_types = FALSE)
    cat("rows:", nrow(d),
        "| clusters:", length(unique(d$cluster)),
        "| comparisons:", paste(unique(d$clone_comparison), collapse = ";"), "\n")
    paths <- c(paths, out)
  }
}

# Exercise the downstream tabulator on the produced paths (+ a deliberate NA) to
# confirm the .drop_missing_paths guard holds.
cat("\n=== tabulate_diffex_clones with an injected NA path ===\n")
res <- tryCatch({
  tabulate_diffex_clones(
    cluster_diffex_clones = as.list(c(paths[1], NA_character_)),
    cluster_xlsx        = "results/_smoke_diffex_per_cluster.xlsx",
    cluster_by_chr_xlsx = "results/_smoke_diffex_per_cluster_by_chr.xlsx",
    total_diffex_clones = list(NA_character_),
    total_xlsx          = "results/_smoke_diffex_total.xlsx",
    total_by_chr_xlsx   = "results/_smoke_diffex_total_by_chr.xlsx")
  "OK (no NA-path crash)"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(res, "\n")
cat("\nSMOKE DONE\n")
