#!/bin/bash
#SBATCH --job-name=rebuild_scna_collages_2p6p
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_scna_collages_2p6p_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_scna_collages_2p6p_%j.err

# Rebuild two_clone_scna_collage_filtered and all_clone_scna_collage_filtered
# over the widened 2p/6p sample scope, then sync to Drive.
#
# WHAT CHANGED (2026-08-27)
#   - config/large_clone_comparisons.yaml regenerated from the current numbat
#     objects. The curated file predated the majority-of-rounds revamp and 0 of
#     its 53 entries still matched (36 mismatch, 10 partial, 4 dead clone id);
#     e.g. SRX10264523 "5_v_4_2p+" was really a chr19 difference, the 2p
#     acquisition being 7_v_4. Backup: large_clone_comparisons_precrevamp_20260827.yaml.bak
#   - new gate scna_collage_samples (NOT paper_retained_samples, which is shared
#     with other target families): 2p 9 samples, 6p 8 samples, 1q/16q unchanged.
#   - new collage_scna_seus_<scna> path targets keyed on that gate.
#   - 6p added to both filtered collage families for the first time.
#
# One tar_make() at a time against _targets_r431 (AGENTS.md) -- check the queue
# before submitting:  squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"

module load r/4.4.1
module load rclone/1.74.3   # needed by hypoxia_rebuilt_gdrive
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Fail before the collage work rather than after it.
command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote configured"; exit 1; }

# install_local, not load_all: collages are built by crew workers in separate R
# processes, which do not inherit a load_all() namespace.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
cat("numbatHelpers installed\n")
' || { echo "FATAL: install failed"; exit 1; }

# Invalidate both filtered families (all four SCNAs incl. the new 6p) so the new
# clone pairs and sample scope actually take effect. These targets carry
# error="null", so a stale cached branch would otherwise persist silently.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
numbatHelpers::init_seu_metadata_db()   # create cell_scores before workers race for it
meta <- tar_meta(fields="name", complete_only=FALSE)
pats <- c("^two_clone_res_collages_filtered_(1q|2p|6p|16q)(_[0-9a-f]+)?$",
          "^all_clone_res_collages_filtered_(1q|2p|6p|16q)(_[0-9a-f]+)?$",
          "^collage_scna_seus_(1q|2p|6p|16q)(_[0-9a-f]+)?$",
          "^scna_collage_samples$",
          "^large_clone_comparisons$")
to_inv <- unique(unlist(lapply(pats, function(p) meta$name[grepl(p, meta$name)])))
if (length(to_inv)) tar_invalidate(any_of(to_inv))
cat("invalidated", length(to_inv), "targets/branches\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

# hypoxia_rebuilt_gdrive is the sync; it also depends on the rest of the hypoxia
# deliverables, which are current after job 11396137.
Rscript run_targets_bg.R \
  two_clone_res_collages_filtered_1q two_clone_res_collages_filtered_2p \
  two_clone_res_collages_filtered_6p two_clone_res_collages_filtered_16q \
  all_clone_res_collages_filtered_1q all_clone_res_collages_filtered_2p \
  all_clone_res_collages_filtered_6p all_clone_res_collages_filtered_16q \
  hypoxia_rebuilt_gdrive

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"rebuild_scna_collages_2p6p ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
