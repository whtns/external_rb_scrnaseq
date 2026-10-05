# Regenerate the Fig. 2 single-sample panels (plot_fig_02) from the ALREADY-BUILT
# t=1e-5 objects, without tar_make().
#
# tar_make(fig_single_sample_panels) would first rebuild every outdated upstream
# target, including the filtered objects that still use the pre-t=1e-5 segment
# letters in config/large_filter_expression.yaml. This reads the built low-hypoxia
# Seurat paths and numbat RDS paths with tar_read() and only redraws the figure.
#
# Writes results/fig_02_<SRX>.pdf for each single_sample_panel_ids sample
# (SRX11133594 is the manuscript's Fig. 2). Reads the held-out samples' outputs;
# writes nothing into output/numbat_*.

suppressPackageStartupMessages({
  library(targets)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")

ids        <- tar_read(single_sample_panel_ids)
seu_paths  <- unlist(tar_read(seus_low_hypoxia))
nb_files   <- tar_read(numbat_rds_files)
simplify   <- tar_read(large_clone_simplifications)

only <- Sys.getenv("FIG02_IDS")
if (nzchar(only)) ids <- strsplit(only, ",")[[1]]

for (id in ids) {
  seu_path <- seu_paths[grepl(id, seu_paths)]
  if (length(seu_path) != 1) { message("skip ", id, ": ", length(seu_path), " seu paths"); next }
  message(Sys.time(), "  ", id, "  <- ", seu_path)
  res <- tryCatch(plot_fig_02(seu_path, nb_files, simplify),
                  error = function(e) { message("  FAILED ", id, ": ", conditionMessage(e)); NA })
  message("  -> ", res)
}
