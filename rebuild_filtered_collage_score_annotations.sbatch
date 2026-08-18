#!/bin/bash
#SBATCH --job-name=rebuild_filtered_collage_scores
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_filtered_collage_scores_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_filtered_collage_scores_%j.err

# Rebuild the two FILTERED collage families with the continuous hypoxia_score /
# mito_score column annotations (see docs/collage_score_annotations.md), plus the
# targets that collect them into figures. The low-hypoxia two-clone family is NOT
# invalidated -- it keeps the cell-cycle-only annotation by design.
#
# One tar_make() at a time against _targets_r431 (AGENTS.md) -- check the queue
# before submitting:  squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"

module load r/4.4.1
module load rclone/1.74.3   # needed by hypoxia_rebuilt_gdrive (bare `rclone` does not resolve here)
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Fail before the 10 h of collage work rather than after it: a missing rclone
# makes hypoxia_rebuilt_gdrive error, and error="null" would cache that as a
# silent NULL.
command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote configured"; exit 1; }

# install_local, not load_all: the collages are built by crew workers in separate
# R processes, which do not inherit a load_all() namespace.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
cat("numbatHelpers installed\n")
' || { echo "FATAL: install failed"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
numbatHelpers::init_seu_metadata_db()   # create cell_scores before the workers race for it
meta <- tar_meta(fields="name", complete_only=FALSE)
pats <- c("^two_clone_res_collages_filtered_(1q|2p|16q)(_[0-9a-f]+)?$",
          "^all_clone_res_collages_filtered_(1q|2p|16q)(_[0-9a-f]+)?$")
to_inv <- unique(unlist(lapply(pats, function(p) meta$name[grepl(p, meta$name)])))
tar_invalidate(any_of(to_inv))
cat("invalidated", length(to_inv), "filtered collage targets/branches\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

# Rebuilds the PDFs in place under results/, then syncs the whole hypoxia_rebuilt
# tree to Drive. hypoxia_rebuilt_gdrive already lists both filtered collage
# families (two_clone_scna_collage_filtered / all_clone_scna_collage_filtered),
# so the annotated PDFs land under their existing subdirectories -- but it also
# depends on the rest of the hypoxia deliverables, which must be up to date or
# they will rebuild here too.
Rscript run_targets_bg.R \
  two_clone_res_collages_filtered_1q two_clone_res_collages_filtered_2p two_clone_res_collages_filtered_16q \
  all_clone_res_collages_filtered_1q all_clone_res_collages_filtered_2p all_clone_res_collages_filtered_16q \
  hypoxia_rebuilt_gdrive

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"rebuild_filtered_collage_scores ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
