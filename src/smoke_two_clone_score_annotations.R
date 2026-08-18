#!/usr/bin/env Rscript
# Smoke test for the hypoxia / mito column annotations on the FILTERED collages:
# one sample, one SCNA, one resolution, both builders. Confirms that
#   1. .ensure_hypoxia_scores() computes hypoxia_score + mito_score on a
#      *_filtered_seu.rds (which carries neither) and CACHES them in cell_scores,
#   2. a second call READS the cache instead of rescoring,
#   3. the two-clone and all-clone filtered collages both render with the
#      annotations drawn as continuous bands in their pinned colours.
suppressPackageStartupMessages({
  source("packages.R")
  library(targets); library(Seurat); library(stringr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")

lcc <- tar_read(large_clone_comparisons)
nbp <- tar_read(numbat_rds_files)
lcs <- tar_read(large_clone_simplifications)
hip <- tar_read(seus_high_hypoxia)

seu_path <- "output/seurat/SRX10264523_filtered_seu.rds"
annots   <- c("hypoxia_score", "mito_score")

# --- 1. cold scoring + cache write ---------------------------------------
init_seu_metadata_db()
con <- numbatHelpers:::connect_hash_db("batch_hashes.sqlite")
DBI::dbExecute(con, "DELETE FROM cell_scores WHERE filepath = ?", params = list(seu_path))
DBI::dbDisconnect(con)
cat("cleared any existing cell_scores rows for", seu_path, "\n")

seu <- readRDS(seu_path)
cat("before:", paste(intersect(c("hypoxia_score", "MT", "mito_score"),
                               colnames(seu@meta.data)), collapse = ", "), "(none expected)\n")
t0 <- Sys.time()
seu <- numbatHelpers:::.ensure_hypoxia_scores(seu, cache_path = seu_path)
cat("cold scoring took", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
for (cl in annots) {
  if (!cl %in% colnames(seu@meta.data)) { cat("!! MISSING after ensure:", cl, "\n"); next }
  v <- seu@meta.data[[cl]]
  cat(sprintf("cold  : %-13s range %.4f .. %.4f (NA %d)\n",
              cl, min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(is.na(v))))
}
cold <- seu@meta.data[, c("hypoxia", "MT", "hypoxia_score")]

# --- 2. warm read from the cache -----------------------------------------
seu2 <- readRDS(seu_path)
t0 <- Sys.time()
seu2 <- numbatHelpers:::.ensure_hypoxia_scores(seu2, cache_path = seu_path)
cat("warm scoring took", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
warm <- seu2@meta.data[, c("hypoxia", "MT", "hypoxia_score")]
cat("cache round-trips identically:",
    isTRUE(all.equal(cold, warm, check.attributes = FALSE)), "\n")
rm(seu, seu2, cold, warm); gc()

# --- 3. both collage builders through their real entry points -------------
two_out <- plot_scna_two_clone_res_collages(
  seu_path, scna_of_interest = "1q", large_clone_comparisons = lcc,
  resolutions = 0.6, nb_paths = nbp, clone_simplifications = lcs,
  score_annotations = annots
)
cat("\ntwo-clone returned:", paste(two_out, collapse = ", "), "\n")
cat("two-clone pdf ok:", length(two_out) == 1 && !is.na(two_out) && file.exists(two_out), "\n")

all_out <- plot_scna_all_clone_res_collages(
  seu_path, scna_of_interest = "1q", resolutions = 0.6,
  high_hypoxia_paths = hip, nb_paths = nbp, clone_simplifications = lcs,
  score_annotations = annots
)
cat("all-clone returned:", paste(all_out, collapse = ", "), "\n")
cat("all-clone pdf ok:", length(all_out) == 1 && !is.na(all_out) && file.exists(all_out), "\n")

cat("\nSMOKE DONE\n")
