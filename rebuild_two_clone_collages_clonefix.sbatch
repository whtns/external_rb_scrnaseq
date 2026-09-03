#!/bin/bash
#SBATCH --job-name=two_clone_clonefix
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/two_clone_clonefix_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/two_clone_clonefix_%j.err

# Re-render the two-clone SCNA collages with the clone-identity fixes in
# numbatHelpers 74ff1fc, rebuild the merged per-SCNA PDFs, and sync to Drive.
#
# WHAT CHANGED (numbatHelpers 74ff1fc)
#   - plot_clone_tree() labelled SUBSET trees with the node's index in the FULL
#     numbat tree instead of its clone id, so the two-clone 6p collage for
#     SRX14116946 (clones 4 and 7) drew a tree reading "8 -> 10" beside panels
#     reading clone 4 / clone 7, and SRX10264523 2p drew "6 -> 13". Subset trees
#     now carry numbat clone ids; whole trees are untouched.
#   - plot_scna_two_clone_res_collages() took the UNION of the clones over every
#     comparison matching the SCNA, so SRX10264526 1q (matching both
#     2_v_1_1q+_16q- and 4_v_3_1q+) built a FOUR-clone collage. Each comparison
#     now gets its own strictly two-clone set, tagged _1q-2v1_ / _1q-4v3_.
#
# EXPECTED OUTPUT CHANGES
#   - every two-clone filtered collage re-renders; only the tree/segment-tree
#     panels and their node colours differ, and only where the pair is not 1/2
#   - SRX10264526 1q: 7 pages replaced by 14 (two comparisons x 7 resolutions)
#     under new names; the old _1q_res*.pdf for that sample are superseded
#   - the four merged per-SCNA PDFs pick the new pages up automatically
#     (cue="always"), so 1q grows from 42 pages to 49
#
# NOT rebuilt here: the hypoxia_low two-clone twins
# (two_clone_res_collages_<scna>) and the all-clone families carry the same tree
# fix but are a separate multi-hour render; they will pick it up when next built.
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

# The collage targets depend on numbatHelpers only through the function body,
# which targets does not hash, so nothing goes stale on its own -- invalidate
# explicitly. These targets carry error="null", so a stale cached branch would
# otherwise persist silently.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
numbatHelpers::init_seu_metadata_db()   # create cell_scores before workers race for it
meta <- tar_meta(fields="name", complete_only=FALSE)
pats <- c("^two_clone_res_collages_filtered_(1q|2p|6p|16q)(_[0-9a-f]+)?$")
to_inv <- unique(unlist(lapply(pats, function(p) meta$name[grepl(p, meta$name)])))
if (length(to_inv)) tar_invalidate(any_of(to_inv))
cat("invalidated", length(to_inv), "targets/branches\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

Rscript run_targets_bg.R \
  two_clone_res_collages_filtered_1q two_clone_res_collages_filtered_2p \
  two_clone_res_collages_filtered_6p two_clone_res_collages_filtered_16q \
  two_clone_res_collages_filtered_merged_1q two_clone_res_collages_filtered_merged_2p \
  two_clone_res_collages_filtered_merged_6p two_clone_res_collages_filtered_merged_16q \
  hypoxia_rebuilt_gdrive

# Read the clone-tree node labels back out of the rendered PDFs: the fix is
# working when these are the numbat clone ids the bar panel names, not 8/10.
module load poppler/23.04.0 2>/dev/null
for f in results/SRX14116946__6p_res0.6_heatmap_phase_scatter_patchwork.pdf \
         results/SRX10264523__2p_res0.6_heatmap_phase_scatter_patchwork.pdf; do
  [ -f "$f" ] && echo "CHECK $f: $(pdftotext -f 1 -l 1 "$f" - 2>/dev/null | grep -cE '^(4|7)$') node labels in {4,7}"
done
ls -1 results/*_1q-*v*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | \
  sed 's/^/SPLIT: /' || echo "SPLIT: none written"
for f in results/two_clone_scna_collage_filtered_merged_*.pdf; do
  [ -f "$f" ] && echo "MERGED: $f ($(stat -c%s "$f") bytes)"
done

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"two_clone_clonefix ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
