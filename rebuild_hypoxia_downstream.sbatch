#!/bin/bash
#SBATCH --job-name=rebuild_hypoxia_downstream
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_downstream_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_downstream_%j.err

# Rebuild the STALE hypoxia-downstream plots/tables (scope b: heatmaps + clone
# trees + diffex + sample_summaries) against the new Jul 2 hypoxia subsetting.
#
# The split logic changed (commit 2f8fbf9) but seus_low/high_hypoxia return
# unchanged path strings, so nothing downstream re-ran (path-string staleness
# trap, see AGENTS.md / memory seu_targets_depend_false_freeze). This job
# force-invalidates the FULL downstream closure of the split, then rebuilds the
# scope-b deliverables. tar_make walks their upstream deps -- all invalidated,
# so every intermediate (hypoxia_seus_1q/2p/6p/16q, integrated_seu_low_hypoxia*,
# numbat_heatmap_plots_low_hypoxia, ...) rebuilds fresh instead of skipping.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Crew workers need the current package. document() first so @export tags added
# since the last NAMESPACE regen take effect (e.g. myreadxl). upgrade="never" so a
# failing dependency upgrade (xfun last run) can't abort the install and leave a
# stale build in which internal helpers like myreadxl are missing.
Rscript -e 'suppressMessages({devtools::document("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet=TRUE); devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE)}); stopifnot(exists("myreadxl", envir=asNamespace("numbatHelpers"), inherits=FALSE)); cat("numbatHelpers installed; myreadxl present\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed -- aborting before tar_make"; exit 1; }

# Invalidate the downstream closure of the split (computed from the DAG, so no
# intermediate is missed). A safety assertion aborts if the closure ever reaches
# a known-upstream target -- guards against an edge-direction mistake nuking
# filtered_seus / cluster_dictionary / the numbat inputs.
Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store = "_targets_r431")
net <- tar_network(targets_only = TRUE, outdated = FALSE, callr_function = NULL)
E <- net$edges  # from = upstream dependency, to = downstream target
seeds <- c("seus_low_hypoxia", "seus_high_hypoxia", "hypoxia_partition_paths")
# base:: qualify all set ops: the targets globals load `conflicted`, which aborts
# on bare union/setdiff/intersect ("found in 4 packages"). That crash previously
# left NOTHING invalidated -> tar_make skipped 918 targets (path-string trap).
ds <- character(0); fr <- seeds
repeat {
  nx  <- unique(E$to[E$from %in% fr])
  new <- base::setdiff(nx, c(ds, seeds))
  if (!length(new)) break
  ds  <- base::union(ds, new); fr <- new
}
forbidden <- c("filtered_seus","unfiltered_seus","cluster_dictionary",
               "cluster_dictionary_per_sample","numbat_rds_files","hypoxia_seus")
bad <- base::intersect(ds, forbidden)
if (length(bad)) stop("SAFETY ABORT: downstream closure reached upstream target(s): ",
                      paste(bad, collapse = ", "))
# Already refreshed by the Jul 2 atomic rebuild / the effect_of_filtering plot job.
fresh <- c("effect_of_filtering","effect_of_filtering_collated","hypoxia_split_log_collated")
to_inv <- base::setdiff(ds, fresh)
tar_invalidate(any_of(to_inv))
cat("Downstream closure size:", length(ds), "| invalidated:", length(to_inv), "\n")
cat(paste(sort(to_inv), collapse = "\n"), "\n")
stopifnot(length(to_inv) > 0)  # a 0-count means the closure walk failed silently
' \
  || { echo "FATAL: invalidation step failed -- aborting before tar_make"; exit 1; }

# Rebuild scope-b deliverables (heatmaps + clone trees + diffex + sample_summaries).
Rscript run_targets_bg.R \
  sample_summaries \
  heatmap_collages_low_hypoxia heatmap_collages_high_hypoxia heatmap_collages_6p \
  numbat_heatmap_plots_low_hypoxia hypoxia_gene_heatmap_low_hypoxia \
  hypoxia_low_integrated_heatmap_collages annotated_heatmap_collages \
  low_hypoxia_clone_tree_files low_hypoxia_clone_trees_segments_files \
  debranched_clone_trees debranched_clone_tree_files \
  all_diffex_clones cis_diffex_clones trans_diffex_clones \
  diffex_6p_g1 fig_1q_cluster_diffex_unintegrated_input fig_16q_cluster_diffex_unintegrated_input
# EXCLUDED targets:
#  - clustree_diffexes (per request): self-contained clustree branch
#    (clustrees -> clustree_tables -> clustree_diffexes -> clustree_*_changes);
#    omitting it also skips the clustree_tables build.
#  - {all,cis,trans}_diffex_clones_for_each_cluster: one sample-branch throws
#    "`x` must be a vector" (degenerate clusters/clone_opt in a fresh low-hypoxia
#    subset), which with error="stop" halted the whole run before the heatmaps/
#    trees/sample_summaries could build. Verified via tar_network that NO other
#    deliverable depends on these, so dropping them unblocks everything else.
#    The per-cluster diffex needs a separate root-cause fix in
#    find_diffex_bw_clones_for_each_cluster (numbatHelpers plot_functions_26.R).
