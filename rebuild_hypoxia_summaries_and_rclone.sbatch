#!/bin/bash
#SBATCH --job-name=rebuild_hypsumm_rclone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypsumm_rclone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypsumm_rclone_%j.err

# Regenerate the phase-labelled hypoxia collages (low/high/whole + annotated),
# recollate hypoxia_summaries from them, then rclone the whole hypoxia deliverable
# set to Google Drive -- all in one tar_make of hypoxia_rebuilt_gdrive (which
# depends on hypoxia_summaries + annotated_heatmap_collages, so the chain cascades).
#
# The collage builders call plot_seu_marker_heatmap (auto_phase_level); a
# numbatHelpers function-body change is not command-tracked by targets, so we
# tar_invalidate the collage + summary targets, then rebuild. rclone module is
# loaded for the deployment="main" rclone target.

module load r/4.4.1 rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
mkdir -p tmp   # low/high collages are tempfiles here that hypoxia_summaries reads

# Ensure crew workers have the current auto-phase code installed.
Rscript -e 'devtools::document("/project2/cobrinik_1090/rpkgs/numbat_helpers")' \
  && Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force = TRUE, upgrade = "never")' \
  || { echo "FATAL: package build/install failed"; exit 1; }
Rscript -e 'stopifnot(exists("assign_auto_phase_clusters", where = asNamespace("numbatHelpers")))' \
  || { echo "FATAL: auto-phase helper not installed"; exit 1; }

# Invalidate the phase-labelled collage targets + the summary so they rebuild
# with the auto-generated cluster names; everything else (clone trees, numbat
# heatmaps, sample_summaries, effect_of_filtering) stays cached.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
targets_re <- "^(heatmap_collages_low_hypoxia|heatmap_collages_high_hypoxia|heatmap_collages_hypoxia|annotated_heatmap_collages|hypoxia_summaries)(_[0-9a-f]+)?$"
meta <- tar_meta(fields = "name", complete_only = FALSE)
to_inv <- meta$name[grepl(targets_re, meta$name)]
if (length(to_inv)) tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "collage/summary targets/branches\n")
' || { echo "FATAL: invalidation errored"; exit 1; }

# tar_make hypoxia_rebuilt_gdrive: rebuilds the invalidated collages -> hypoxia_summaries,
# then (deployment="main") rclones the hypoxia deliverable set to Google Drive.
Rscript run_targets_bg.R hypoxia_rebuilt_gdrive
