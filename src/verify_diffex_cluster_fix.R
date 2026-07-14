# Verify the find_diffex_bw_clones_for_each_cluster() fix on the objects that
# actually broke: call it exactly as the fig_16q / fig_1q targets do and confirm it
# now returns a real diffex CSV path instead of throwing "`x` must be a vector".
options(warn = 1)
suppressMessages({
  library(targets)
  library(tidyverse)
  library(Seurat)
  library(seuratTools)
})
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
tar_config_set(store = "_targets_r431")

numbat_rds_files       <- tar_read_raw("numbat_rds_files")
# read the yaml directly: the target was invalidated when it was re-keyed, and this
# is the artifact under test anyway
large_clone_comparisons <- yaml::read_yaml("config/large_clone_comparisons.yaml")
cluster_dictionary     <- tar_read_raw("cluster_dictionary")
cat("clone-comparison keys:", length(large_clone_comparisons),
    "| all SRX:", all(startsWith(names(large_clone_comparisons), "SRX")), "\n")

for (spec in list(list(tgt = "hypoxia_seus_16q", scna = "16q-"),
                  list(tgt = "hypoxia_seus_1q",  scna = "1q+"))) {
  paths <- unlist(tar_read_raw(spec$tgt))
  cat("\n=====", spec$tgt, "scna:", spec$scna, "(", length(paths), ")\n")
  for (f in paths) {
    out <- tryCatch(
      find_diffex_bw_clones_for_each_cluster(
        f, numbat_rds_files, large_clone_comparisons, cluster_dictionary,
        location = "all", scna_of_interest = spec$scna),
      error = function(e) paste0("!! STILL ERRORS: ", conditionMessage(e)))
    n <- if (!is.na(out) && !grepl("^!!", out) && file.exists(out)) {
      nrow(readr::read_csv(out, show_col_types = FALSE))
    } else NA_integer_
    cat(sprintf("  %-42s -> %s  (rows: %s)\n", basename(f),
                if (is.na(out)) "NA (skipped)" else out, n))
  }
}
cat("\nVERIFY DONE\n")
