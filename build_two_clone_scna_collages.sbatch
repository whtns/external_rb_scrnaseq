#!/bin/bash
#SBATCH --job-name=two_clone_scna_collages
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/two_clone_scna_collages_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/two_clone_scna_collages_%j.err

# Rebuild the two-clone (acquiring vs immediately-preceding) SCNA collages for
# 1q+/2p+/16q- across paper_retained samples, as a TRUE ZOOM of the full-population
# low-hypoxia collages.
#
# plot_scna_two_clone_res_collages() reads the SCT_snn_res.* columns the hypoxia
# split persisted on each low object (PCA recomputed on the low-hypoxia population
# upstream) -- NO PCA / NO clustering here. It then hands the FULL object to
# plot_seu_marker_heatmap(display_cells=<the two clones' cells>): cluster order,
# phase labels, and marker-gene rows are all derived on the full object and only
# then are the panels restricted to the two clones. So each two-clone panel carries
# EXACTLY the same cluster IDs, marker rows, and order as its full-population twin.
# Capped at res 0.2..0.8 to mirror the full low-hypoxia sweep, so every two-clone
# panel has a same-resolution full-population counterpart. Output:
#   results/<SRX>_hypoxia_low_seu.rds__<scna>_res<R>_heatmap_phase_scatter_patchwork.pdf
# Re-synced to Drive via hypoxia_rebuilt_gdrive under
#   two_clone_scna_collage/{paper_retained,other}/.
# (6p omitted: no paper_retained sample carries a 6p+ subclone.)

module load r/4.4.1
module load poppler/23.04.0
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers. NAMESPACE is HAND-MAINTAINED; verify the exports on
# the EXPORTED surface. Also verify plot_seu_marker_heatmap gained `display_cells`
# and the two-clone fn no longer takes `recompute_pca` (guard against stale install).
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
ns <- asNamespace("numbatHelpers")
ex <- getNamespaceExports("numbatHelpers")
need <- c("plot_scna_two_clone_res_collages","plot_seu_marker_heatmap")
miss <- setdiff(need, ex); if (length(miss)) stop("NOT EXPORTED: ", paste(miss, collapse=", "))
if (!"display_cells" %in% names(formals(get("plot_seu_marker_heatmap", ns))))
  stop("STALE install: plot_seu_marker_heatmap lacks display_cells")
if ("recompute_pca" %in% names(formals(get("plot_scna_two_clone_res_collages", ns))))
  stop("STALE install: plot_scna_two_clone_res_collages still has recompute_pca")
cat("numbatHelpers installed; display_cells present; recompute_pca gone\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed -- aborting before tar_make"; exit 1; }

# Cap is now 0.8, so remove the res1.0 two-clone PDFs from disk BEFORE the build
# (rclone copy never deletes, and a 0.8-capped rebuild won't touch a res1.0 file).
echo "=== deleting res1.0 two-clone PDFs from disk ==="
ls results/*_hypoxia_low_seu.rds__*_res1_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -E "__(1q|2p|16q)_res1_" | wc -l
rm -f results/*_hypoxia_low_seu.rds__1q_res1_heatmap_phase_scatter_patchwork.pdf \
      results/*_hypoxia_low_seu.rds__2p_res1_heatmap_phase_scatter_patchwork.pdf \
      results/*_hypoxia_low_seu.rds__16q_res1_heatmap_phase_scatter_patchwork.pdf

# The two-clone targets call package-internal functions (plot_scna_two_clone_res_collages
# -> plot_seu_marker_heatmap), which targets does NOT hash for invalidation, and the
# target COMMAND is unchanged this round -- so invalidate the branches explicitly or
# tar_make would skip them. hypoxia_rebuilt_gdrive is cue="always".
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("two_clone_res_collages_1q","two_clone_res_collages_2p","two_clone_res_collages_16q")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "two-clone targets/branches\n")
stopifnot(length(to_inv) > 0)
' || { echo "FATAL: invalidation failed -- aborting"; exit 1; }

Rscript run_targets_bg.R \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  hypoxia_rebuilt_gdrive

# rclone copy never deletes, so purge the res1.0 two-clone PDFs from Drive too
# (restricted to the two_clone_scna_collage tree so nothing else is touched).
echo "=== deleting res1.0 two-clone PDFs from Drive ==="
rclone delete gdrive:rb_scrnaseq/hypoxia_rebuilt/two_clone_scna_collage \
  --include "*_res1_heatmap_phase_scatter_patchwork.pdf" --verbose 2>&1 | tail -5

# Slack ping on completion (best-effort).
if [ -f "$HOME/.slack_webhook" ]; then
  N=$(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -E "__(1q|2p|16q)_res" | wc -l)
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"two_clone_scna_collages ($SLURM_JOB_ID) done: $N two-clone SCNA collages (lock-to-full zoom) on disk, synced to Drive\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
