#!/bin/bash
#SBATCH --job-name=hyp_downstream
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_downstream_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_downstream_%j.err

# Rebuilds the split-dependent low-hypoxia deliverables with the new phase gate,
# then rclones them to Google Drive. These targets key off deterministic path
# STRINGS from seus_low_hypoxia (not file content), so the split change does NOT
# re-key them (the seu_targets_depend_false_freeze issue) -- tar_invalidate below
# forces the recompute. Shared intermediates (cluster_orders,
# large_clone_simplifications*, numbat plots) are file/numbat-derived, not from
# the hypoxia subset, so they are intentionally NOT invalidated.
#
# Scope (per request): clone trees, numbat heatmaps, sample_summaries,
# effect_of_filtering, annotated collages -> hypoxia_rebuilt_gdrive.

module load r/4.4.1 rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

INVAL='c("low_hypoxia_clone_tree_files","low_hypoxia_clone_trees_segments_files","numbat_heatmap_plots_low_hypoxia","annotated_heatmap_collages","effect_of_filtering","sample_summaries","hypoxia_rebuilt_gdrive")'
Rscript -e "library(targets); tar_config_set(store='_targets_r431'); tar_invalidate(any_of(${INVAL})); cat('invalidated\n')" || { echo "FATAL: invalidation errored"; exit 1; }

Rscript run_targets_bg.R hypoxia_rebuilt_gdrive
