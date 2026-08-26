#!/bin/bash
#SBATCH --job-name=cluster_diffex_figs
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/cluster_diffex_figs_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/cluster_diffex_figs_%j.err

# Build the per-cluster clone diffex figures (Fig. 7b / 8b) now that BOTH stacked
# bugs are fixed:
#   numbatHelpers cc98d4e  - find_diffex_bw_clones_for_each_cluster() derives the
#                            missing `clusters` column instead of throwing
#                            "`x` must be a vector" on every low-hypoxia input
#   project e230a98        - large_clone_comparisons.yaml re-keyed SRR -> SRX, so
#                            the lookup stops silently returning NULL
#
# Verified offline (job 10217640): 8 of 9 inputs yield real diffex CSVs, 61k-102k
# rows. These figures have never actually built with data.
#
# large_clone_comparisons was invalidated when the yaml was re-keyed, so the whole
# diffex chain rebuilds from it.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Crew workers load the INSTALLED package. Two DIFFERENT checks, because the two
# kinds of function have different requirements:
#   - functions the PIPELINE calls must be EXPORTED (it reaches them via
#     library(numbatHelpers)). exists(..., asNamespace=) is not sufficient: it also
#     sees unexported internals, so it passes even when the pipeline cannot call
#     them -- that false pass burned a whole run earlier today.
#   - functions called only from INSIDE the package (assign_auto_phase_clusters,
#     used by find_diffex_bw_clones_for_each_cluster) need only exist in the
#     namespace; package internals are visible to their siblings. Demanding an
#     export for these is too strict and wrongly aborts a valid build.
Rscript -e 'suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE)); ns <- asNamespace("numbatHelpers"); ex <- getNamespaceExports("numbatHelpers"); need_export <- c("find_diffex_bw_clones_for_each_cluster"); need_exist <- c("assign_auto_phase_clusters","auto_phase_level"); m1 <- setdiff(need_export, ex); m2 <- Filter(function(f) !exists(f, envir=ns, inherits=FALSE), need_exist); if (length(m1)) stop("NOT EXPORTED: ", paste(m1, collapse=", ")); if (length(m2)) stop("MISSING FROM NAMESPACE: ", paste(m2, collapse=", ")); cat("numbatHelpers installed; exports + internals present\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed"; exit 1; }

Rscript run_targets_bg.R \
  fig_1q_cluster_diffex_unintegrated fig_16q_cluster_diffex_unintegrated \
  fig_1q_cluster_diffex_integrated fig_16q_cluster_diffex_integrated
