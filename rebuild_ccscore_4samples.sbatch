#!/bin/bash
#SBATCH --job-name=rebuild_ccscore_4samples
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_ccscore_4samples_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_ccscore_4samples_%j.err

# Regenerate cell-cycle scores for the 4 flagged samples by re-running the
# hypoxia chain (hypoxia_seus -> hypoxia_partition_paths -> seus_low/high ->
# collages/summaries) for ONLY their branches, with the numbatHelpers fix that
# scores CellCycleScoring on the gene assay instead of SCT.
#
# hypoxia_seus / hypoxia_partition_paths use cue(command=FALSE), so a function
# body change does NOT auto-invalidate; we tar_invalidate the 4 samples'
# branches explicitly, then dependency-tracking cascades to seus_low/high and
# the plots. Crew heavy workers need the installed package, so install_local.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Rebuild + install numbatHelpers so crew workers pick up the gene-assay fix.
Rscript -e 'devtools::document("/project2/cobrinik_1090/rpkgs/numbat_helpers")' \
  && Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force = TRUE, upgrade = "never")' \
  || { echo "FATAL: package build/install failed"; exit 1; }

# Verify the fix is actually in the installed package before proceeding.
Rscript -e 'src <- deparse(numbatHelpers::load_and_save_hypoxia_score); if (!any(grepl("score_assay", src))) { cat("FATAL: gene-assay cc fix not in installed package\n"); quit(status = 1) }' \
  || exit 1

# Invalidate ONLY the 4 samples' branches across the hypoxia chain, then tar_make.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
samples <- c("SRX10831282","SRX10831286","SRX11133589","SRX22868105")
pat <- paste(samples, collapse = "|")

# target -> the mapping variable it branches over (same per-sample order)
chain <- c(
  hypoxia_seus                          = "filtered_seus",
  hypoxia_threshold_per_sample          = "hypoxia_seus",
  hypoxia_score_plots                   = "hypoxia_seus",
  hypoxia_partition_paths               = "hypoxia_seus",
  seus_low_hypoxia                      = "hypoxia_partition_paths",
  seus_high_hypoxia                     = "hypoxia_partition_paths",
  heatmap_collages_hypoxia              = "hypoxia_seus",
  numbat_heatmap_plots_low_hypoxia      = "seus_low_hypoxia",
  annotated_heatmap_collages            = "seus_low_hypoxia",
  heatmap_collages_low_hypoxia          = "seus_low_hypoxia",
  heatmap_collages_high_hypoxia         = "seus_high_hypoxia",
  hypoxia_gene_heatmap_low_hypoxia      = "seus_low_hypoxia",
  low_hypoxia_clone_tree_files          = "seus_low_hypoxia",
  low_hypoxia_clone_trees_segments_files= "seus_low_hypoxia",
  hypoxia_summaries                     = "hypoxia_seus"
)

for (i in seq_along(chain)) {
  tgt <- names(chain)[i]; mv <- chain[[i]]
  vals <- tryCatch(tar_read_raw(mv), error = function(e) NULL)
  if (is.null(vals)) { message("skip ", tgt, " (cannot read ", mv, ")"); next }
  if (!is.list(vals)) vals <- as.list(vals)
  idx <- which(vapply(vals, function(v) any(grepl(pat, unlist(v))), logical(1)))
  if (length(idx) == 0) { message("skip ", tgt, " (no matching branches in ", mv, ")"); next }
  bn <- tryCatch(tar_branch_names_raw(tgt, idx), error = function(e) character(0))
  bn <- bn[!is.na(bn)]
  if (length(bn) == 0) { message("skip ", tgt, " (no branch names)"); next }
  tar_invalidate(any_of(bn))
  message(sprintf("%-38s invalidated %d branches (idx %s)", tgt, length(bn), paste(idx, collapse=",")))
}
' || { echo "FATAL: invalidation step errored"; exit 1; }

# Build the cell-cycle-affected deliverables for the 4 samples. Dependency
# tracking pulls the invalidated hypoxia_seus/partition/subset branches through.
Rscript run_targets_bg.R \
  hypoxia_score_plots \
  numbat_heatmap_plots_low_hypoxia \
  annotated_heatmap_collages \
  heatmap_collages_low_hypoxia \
  heatmap_collages_high_hypoxia \
  hypoxia_gene_heatmap_low_hypoxia \
  low_hypoxia_clone_tree_files \
  low_hypoxia_clone_trees_segments_files \
  hypoxia_summaries
