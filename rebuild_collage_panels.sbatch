#!/bin/bash
#SBATCH --job-name=rebuild_collage_panels
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_collage_panels_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_collage_panels_%j.err

# Rebuild every collage-family target touched by issues #34/#35/#37/#38/#36:
#   #34 stacked-bar fill colours now match the clone colours
#   #35 heatmap column-split labels drawn diagonally on top
#   #37 segment-tree panel restored beside the clone tree
#   #38 clustree panel added
#   #36 two-clone collages ALSO built on the standard filtered objects
# and re-sync the hypoxia deliverables to Drive.
#
# NEVER run two tar_make() jobs against _targets_r431 at once -- check the queue
# for run_targets_bg.R / crew-worker jobs before submitting (AGENTS.md).

module load r/4.4.1
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers. NAMESPACE is hand-maintained; this change adds no
# new EXPORT (the new helpers are internal ::: ), only new args on already-exported
# functions -- verify those args are present against the installed copy.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
fm <- names(formals(numbatHelpers::plot_seu_marker_heatmap))
need <- c("segment_tree","clustree","column_label_rot")
miss <- setdiff(need, fm)
if (length(miss)) stop("plot_seu_marker_heatmap lacks args: ", paste(miss, collapse=", "), " -- stale install")
cat("numbatHelpers installed; collage-panel args present\n")
' || { echo "FATAL: numbatHelpers install/verify failed"; exit 1; }

# Invalidate the collage targets + their dynamic branches so every sample re-runs.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
pats <- c("^heatmap_collages_low_hypoxia(_by_res)?(_[0-9a-f]+)?$",
          "^heatmap_collages_hypoxia(_[0-9a-f]+)?$",
          "^annotated_heatmap_collages(_[0-9a-f]+)?$",
          "^two_clone_res_collages_(1q|2p|16q)(_[0-9a-f]+)?$",
          "^two_clone_res_collages_filtered_(1q|2p|16q)(_[0-9a-f]+)?$",
          "^filtered_scna_seus_(1q|2p|6p|16q)(_[0-9a-f]+)?$",
          "^collages_(2p|6p|16q)(_[0-9a-f]+)?$")
to_inv <- unique(unlist(lapply(pats, function(p) meta$name[grepl(p, meta$name)])))
tar_invalidate(any_of(to_inv))
cat("invalidated", length(to_inv), "collage targets/branches\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

# Build the collages + the Drive sync (hypoxia_rebuilt_gdrive stages the two-clone
# filtered set under two_clone_scna_collage_filtered/).
Rscript run_targets_bg.R \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  two_clone_res_collages_filtered_1q two_clone_res_collages_filtered_2p two_clone_res_collages_filtered_16q \
  heatmap_collages_low_hypoxia_by_res hypoxia_rebuilt_gdrive

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"rebuild_collage_panels ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
