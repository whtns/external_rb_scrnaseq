# Smoke test for the collage-panel work behind issues #34/#35/#37/#38:
#   #34 stacked-bar fill colours match the clone colours in the UMAP / heatmap
#   #35 heatmap column-split labels drawn diagonally on top
#   #37 segment-tree panel restored alongside the clone tree
#   #38 clustree panel added
# Rebuilds one two-clone collage (SRX10264523, 1q, res 0.4) with the tree panels
# on, then reports which panels landed and whether the bar colours agree with the
# clone palette.
suppressPackageStartupMessages({
  source("packages.R")
  library(targets); library(Seurat); library(dplyr)
  devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
})
tar_config_set(store = "_targets_r431")

p    <- "output/seurat/SRX10264523_hypoxia_low_seu.rds"
lcc  <- tar_read(large_clone_comparisons)
lcs  <- tar_read(large_clone_simplifications)
nbs  <- tar_read(numbat_rds_files)

cat("===== .bar_fill_from_clone() unit check =====\n")
# clone 2 acquired 1q, clone 1 precedes it. scna_status levels put the acquiring
# clone FIRST, the clone factor sorts 1 then 2 -- the exact ordering mismatch
# that flipped the colours.
md <- data.frame(
  clone = factor(c("1", "1", "2", "2")),
  scna_status = factor(c("preceding (clone 1)", "preceding (clone 1)",
                         "1q+ (clone 2)", "1q+ (clone 2)"),
                       levels = c("1q+ (clone 2)", "preceding (clone 1)")))
got  <- numbatHelpers:::.bar_fill_from_clone(md, "scna_status")
want <- stats::setNames(scales::hue_pal()(2), c("1", "2"))
cat("  clone palette      :", paste(names(want), want, sep = "=", collapse = "  "), "\n")
cat("  bar fill_colors    :", paste(names(got), got, sep = "=", collapse = "  "), "\n")
ok <- identical(base::unname(got[["preceding (clone 1)"]]), base::unname(want[["1"]])) &&
      identical(base::unname(got[["1q+ (clone 2)"]]),      base::unname(want[["2"]]))
cat("  clone 1 -> its own colour, clone 2 -> its own colour:",
    if (ok) "OK" else "!! STILL FLIPPED", "\n")
cat("  (default scale would have given the acquiring clone 2 the colour",
    want[["1"]], "-- clone 1's)\n")

cat("\n===== rebuild two-clone collage with all panels =====\n")
out <- plot_scna_two_clone_res_collages(
  p,
  scna_of_interest        = "1q",
  large_clone_comparisons = lcc,
  resolutions             = 0.4,
  nb_paths                = nbs,
  clone_simplifications   = lcs
)
for (f in out) {
  if (is.na(f)) { cat("  -> NA\n"); next }
  sz <- if (file.exists(f)) file.info(f)$size else NA
  cat(sprintf("  %s | %s bytes | %s\n", basename(f), format(sz, big.mark = ","),
              if (!is.na(sz) && sz > 20000) "OK" else "!! TOO SMALL"))
}

cat("\nSMOKE DONE\n")
