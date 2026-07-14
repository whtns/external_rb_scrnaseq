# Root-cause the "`x` must be a vector" failure in
# find_diffex_bw_clones_for_each_cluster() (numbatHelpers plot_functions_26.R:38).
#
# Hypothesis: line 38-39
#   myclusters <- sort(unique(seu@meta.data[["clusters"]])) %>% set_names(.)
# When the object has no "clusters" metadata column, seu@meta.data[["clusters"]]
# is NULL -> sort(unique(NULL)) is NULL -> rlang::set_names(NULL) throws
# "`x` must be a vector". Same footgun as the nb_paths set_names(NULL) already
# fixed in plot_functions_19.R.
#
# Check the actual inputs to the failing targets (hypoxia_seus_16q / _1q).
suppressMessages({
  library(targets)
  library(Seurat)
})
tar_config_set(store = "_targets_r431")

for (tgt in c("hypoxia_seus_16q", "hypoxia_seus_1q")) {
  paths <- tryCatch(unlist(tar_read_raw(tgt)), error = function(e) NULL)
  cat("\n=====", tgt, "(", length(paths), "objects )\n")
  for (f in paths) {
    if (!file.exists(f)) { cat("  MISSING ", f, "\n"); next }
    s  <- readRDS(f)
    md <- s@meta.data
    cl <- md[["clusters"]]
    cat(sprintf(
      "  %-42s cells=%-6d clusters_col=%-5s n_clusters=%-3s clone_opt=%-5s\n",
      basename(f), ncol(s), !is.null(cl),
      if (is.null(cl)) "NA" else length(unique(cl)),
      "clone_opt" %in% colnames(md)))
    if (is.null(cl)) {
      cat("    -> NO 'clusters' column: sort(unique(NULL)) = NULL,",
          "set_names(NULL) will throw \"`x` must be a vector\"\n")
      cat("    resolution-like cols: ",
          paste(grep("_res\\.|clusters", colnames(md), value = TRUE),
                collapse = ", "), "\n")
    }
    rm(s); gc(verbose = FALSE)
  }
}
cat("\ndone\n")
