#!/bin/bash
#SBATCH --job-name=diag_cell_depth
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_cell_depth_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diag_cell_depth_%j.err

# github #43: the depth diagnostic, plus a render check on the new panel.
#
# Both read batch_hashes.sqlite and write only into results/. No Seurat object
# is loaded and no targets store is touched, so this is cheap and safe to run
# alongside nothing in particular.
#
# load_all(), not install_local(): nothing here runs on a crew worker.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

echo "===================== diagnostic ====================="
Rscript src/diag_cell_depth_provenance.R || exit 1

echo
echo "================== depth panel render =================="
Rscript -e '
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
# the low-depth outlier, a mid sample, and a high one
for (s in c("SRX10031191", "SRX10264519", "SRX11133587")) {
  p <- plot_cell_depth_panel(s, out_dir = "results/cell_depth_SMOKE")
  d <- get_sample_cell_depth(s)
  cat(sprintf("%s: %s  (n=%d, median=%.0f)\n", s,
              if (is.na(p)) "FAILED" else basename(p),
              nrow(d), if (nrow(d)) median(d$nCount_gene) else NA))
}'
