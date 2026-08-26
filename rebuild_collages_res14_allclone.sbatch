#!/bin/bash
#SBATCH --job-name=rebuild_res14_allclone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_res14_allclone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_res14_allclone_%j.err

# Rebuild collages for four changes made this session:
#   1. Fisher caption on the clone-distribution bars now str_wrap()s (plot_functions_2)
#   2. two-clone sweep extended 0.8 -> 1.4 (pipeline_targets_integration)
#   3. clone-tree / segment-tree / clustree panels shrunk 2 bands -> 1 (plot_functions_19)
#   4. NEW all-clone filtered collages, hypoxia-removed cells marked on the clustree
#      (node colour = removed fraction) and a removed-cell UMAP band
#      (plot_scna_all_clone_res_collages + all_clone_res_collages_filtered_* targets)
#
# Changes 1 and 3 live inside numbatHelpers functions, which targets does NOT track,
# so the collage families that render through plot_seu_marker_heatmap must be
# invalidated to pick them up (change 2 alters the target command and would
# invalidate on its own; the all-clone targets are new and build fresh).
#
# NEVER run two tar_make() jobs against _targets_r431 at once -- check the queue
# for run_targets_bg.R / crew-worker jobs before submitting (AGENTS.md).

module load r/4.4.1
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers and verify the new/threaded args are present, so a
# stale install fails loudly instead of silently rendering the old panels.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
fm <- names(formals(numbatHelpers::plot_seu_marker_heatmap))
need <- c("segment_tree","clustree","column_label_rot","mark_removed_col")
miss <- setdiff(need, fm)
if (length(miss)) stop("plot_seu_marker_heatmap lacks args: ", paste(miss, collapse=", "), " -- stale install")
if (!exists("plot_scna_all_clone_res_collages", where=asNamespace("numbatHelpers")))
  stop("plot_scna_all_clone_res_collages missing -- stale install")
cat("numbatHelpers installed; all-clone + marking args present\n")
' || { echo "FATAL: numbatHelpers install/verify failed"; exit 1; }

# Invalidate every collage family that renders through plot_seu_marker_heatmap, so
# the smaller-panel + wrapped-caption changes propagate everywhere they appear
# (same set as rebuild_collage_panels.sbatch). filtered_scna_seus_* is included so
# both the two-clone and the new all-clone filtered collages re-subset cleanly.
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

# Build the two-clone (res 1.4) + NEW all-clone collages and re-sync to Drive.
# hypoxia_rebuilt_gdrive depends on all of them, so it both builds and stages them.
Rscript run_targets_bg.R \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  two_clone_res_collages_filtered_1q two_clone_res_collages_filtered_2p two_clone_res_collages_filtered_16q \
  all_clone_res_collages_filtered_1q all_clone_res_collages_filtered_2p all_clone_res_collages_filtered_16q \
  heatmap_collages_low_hypoxia_by_res hypoxia_rebuilt_gdrive
