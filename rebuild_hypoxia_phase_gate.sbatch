#!/bin/bash
#SBATCH --job-name=hyp_phase_gate
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_phase_gate_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_phase_gate_%j.err

# Reruns the cluster-based hypoxia split with the direct phase-purity gate
# (max_dom_frac = 0.70) + marker gate (min_hyp_markers = 2) at resolutions
# 0.2/0.6/1.0, then rebuilds the hypoxia_summaries review deliverable.
#
# hypoxia_partition_paths uses cue(command = FALSE), so the code change does NOT
# auto-invalidate it; tar_invalidate below forces the recompute. seus_low/high
# and the split-dependent collages are invalidated so they pick up the new
# subsets.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

TARGETS_TO_INVALIDATE='c("hypoxia_partition_paths","seus_low_hypoxia","seus_high_hypoxia","heatmap_collages_low_hypoxia","heatmap_collages_high_hypoxia","hypoxia_gene_heatmap_low_hypoxia","silhouette_plots","hypoxia_split_log_collated","hypoxia_summaries")'
Rscript -e "library(targets); tar_config_set(store = '_targets_r431'); tar_invalidate(any_of(${TARGETS_TO_INVALIDATE})); cat('invalidated\n')" && \
Rscript run_targets_bg.R hypoxia_summaries hypoxia_split_log_collated
