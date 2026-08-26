#!/bin/bash
#SBATCH --job-name=figs_and_tables
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/figures_and_tables_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/figures_and_tables_%j.err

# Full figures_and_tables build, with the package installed FIRST.
#
# build_figures_and_tables.sbatch does NOT install numbatHelpers, so it runs the
# crew workers against whatever was last installed. Job 10226497 started that way
# and would have rebuilt Fig. 7b with the stale plot_fig_07_08() -- i.e. a NULL
# figure from a green pipeline. Install + verify before tar_make.
#
# Fixes now in play:
#   plot_functions_26  derive the missing `clusters` column (was: set_names(NULL))
#   plot_functions_14  plot_fig_07_08() drops NA/missing CSV paths. A sample with
#                      no clone comparison for the SCNA (SRX10831287 has only 11a+,
#                      no 1q+) legitimately returns NA_character_; read_csv() was
#                      taking that NA as a filename and erroring the whole target.
#   config yaml        large_clone_comparisons re-keyed SRR -> SRX
#
# Fig. 8b (fig_16q_cluster_diffex_unintegrated.pdf) already built for the first time
# ever in job 10226488. Fig. 7b (1q) is what the plot_fig_07_08 fix unblocks.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e 'suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE)); ex <- getNamespaceExports("numbatHelpers"); need <- c("plot_fig_07_08","find_diffex_bw_clones_for_each_cluster"); miss <- setdiff(need, ex); if (length(miss)) stop("NOT EXPORTED: ", paste(miss, collapse=", ")); cat("numbatHelpers installed; exports present\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed"; exit 1; }

Rscript run_targets_bg.R figures_and_tables
