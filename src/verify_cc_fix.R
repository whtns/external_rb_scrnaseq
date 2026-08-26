# Verify the #36 CC-score fix: SRX10264523 filtered_seu previously failed the
# two-clone collage with "'G2M.Score','S.Score' not found". Render one res and
# confirm a PDF lands.
suppressPackageStartupMessages({
  source("packages.R"); library(targets); library(Seurat)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")
lcc <- tar_read(large_clone_comparisons); lcs <- tar_read(large_clone_simplifications)
nbs <- tar_read(numbat_rds_files)
out <- plot_scna_two_clone_res_collages(
  "output/seurat/SRX10264523_filtered_seu.rds",
  scna_of_interest = "1q", large_clone_comparisons = lcc,
  resolutions = 0.4, nb_paths = nbs, clone_simplifications = lcs)
cat("returned:", paste(out, collapse=", "), "\n")
ok <- length(out)==1 && !is.na(out) && file.exists(out) && file.info(out)$size > 20000
cat(if (ok) "CC_FIX_OK\n" else "CC_FIX_FAIL\n")
