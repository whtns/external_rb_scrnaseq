#!/bin/bash
#SBATCH --job-name=rebuild_two_clone_collect
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_two_clone_collect_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_two_clone_collect_%j.err

# Re-run ONLY the two-clone collage families + the Drive sync, so the filtered
# two-clone collages (#36) are actually COLLECTED. They already render to disk,
# but the pre-fix target reconstructed the wrong PDF name and returned NA. The
# numbatHelpers fix (fe00b1f) makes plot_seu_marker_heatmap return the real path;
# this rebuild picks that up. Other collage families are already correct and are
# NOT invalidated here.
#
# One tar_make() at a time against _targets_r431 (AGENTS.md) -- check the queue
# before submitting.

module load r/4.4.1
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
cat("numbatHelpers installed\n")
' || { echo "FATAL: install failed"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
pats <- c("^two_clone_res_collages_(1q|2p|16q)(_[0-9a-f]+)?$",
          "^two_clone_res_collages_filtered_(1q|2p|16q)(_[0-9a-f]+)?$")
to_inv <- unique(unlist(lapply(pats, function(p) meta$name[grepl(p, meta$name)])))
tar_invalidate(any_of(to_inv))
cat("invalidated", length(to_inv), "two-clone targets/branches\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

Rscript run_targets_bg.R \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  two_clone_res_collages_filtered_1q two_clone_res_collages_filtered_2p two_clone_res_collages_filtered_16q \
  hypoxia_rebuilt_gdrive

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"rebuild_two_clone_collect ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
