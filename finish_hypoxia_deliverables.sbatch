#!/bin/bash
#SBATCH --job-name=finish_hypoxia_deliv
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/finish_hypoxia_deliverables_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/finish_hypoxia_deliverables_%j.err

# Finish the two deliverables that rebuild_hypoxia_downstream (job 10207002) did
# not reach. That run rebuilt the new lowest-stable-resolution hypoxia split for
# all 33 samples and every heatmap/clone-tree/collage deliverable, then died on a
# branch of fig_16q_cluster_diffex_unintegrated_input:
#
#   Error: `x` must be a vector
#
# That target calls find_diffex_bw_clones_for_each_cluster() (plot_functions_26.R),
# the SAME known-broken per-cluster diffex that rebuild_hypoxia_downstream.sbatch
# already excludes {all,cis,trans}_diffex_clones_for_each_cluster for -- it throws
# on degenerate clusters/clone_opt in a fresh low-hypoxia subset. fig_16q (and
# fig_1q) were left in the list and, having error="stop", halted the pipeline
# before sample_summaries / hypoxia_split_log_collated could build.
#
# This is a PRE-EXISTING bug, unrelated to the split change, and needs a root-cause
# fix in find_diffex_bw_clones_for_each_cluster. Until then, build the remaining
# deliverables without the fig_1q/fig_16q per-cluster diffex figures.
#
# NO invalidation here: everything upstream is already built and cached by 10207002.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript run_targets_bg.R sample_summaries hypoxia_split_log_collated
