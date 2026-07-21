#!/bin/bash
#SBATCH --job-name=two_clone_and_trees
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/two_clone_scna_collages_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/two_clone_scna_collages_%j.err

# Rebuild the two-clone SCNA collages (1q/2p/16q, paper_retained, 0.2..0.8, true
# ZOOM of the full-population low-hypoxia collages) AND refresh the low-hypoxia
# clone + segment trees.
#
# 1) two_clone collages: plot_scna_two_clone_res_collages() reads the persisted
#    SCT_snn_res.* columns (no PCA/clustering) and hands the FULL object to
#    plot_seu_marker_heatmap(display_cells=<the two clones>, bar_var="scna_status"),
#    so the panels are a true zoom AND the stacked-bar panel is now labeled by the
#    SCNA of interest -- "1q+ (clone 2)" vs "preceding (clone 1)" -- instead of raw
#    numbat clone ids.
# 2) low-hypoxia trees: low_hypoxia_clone_tree_files /
#    low_hypoxia_clone_trees_segments_files map over seus_low_hypoxia, which returns
#    path STRINGS -- so the 2026-07-17 low-object rebuild did NOT invalidate them
#    (they were last built 2026-07-13). Invalidate explicitly so they rebuild from
#    the current retained-cell barcodes. Output: results/clone_trees/*.pdf

module load r/4.4.1
module load poppler/23.04.0
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers. NAMESPACE is HAND-MAINTAINED; verify exports and
# that plot_seu_marker_heatmap gained display_cells + bar_var (guard stale install).
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
ns <- asNamespace("numbatHelpers")
ex <- getNamespaceExports("numbatHelpers")
need <- c("plot_scna_two_clone_res_collages","plot_seu_marker_heatmap","save_clone_tree_from_path")
miss <- setdiff(need, ex); if (length(miss)) stop("NOT EXPORTED: ", paste(miss, collapse=", "))
fm <- names(formals(get("plot_seu_marker_heatmap", ns)))
for (a in c("display_cells","bar_var")) if (!a %in% fm) stop("STALE install: plot_seu_marker_heatmap lacks ", a)
cat("numbatHelpers installed; display_cells + bar_var present; exports ok\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed -- aborting before tar_make"; exit 1; }

# targets does NOT hash package-internal functions, and the two-clone target
# command is unchanged -- so invalidate the two-clone targets explicitly. The
# low-hypoxia tree targets map over seus_low_hypoxia (path strings), so the
# 07-17 low-object rebuild left them "up to date" too -- invalidate them as well.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("two_clone_res_collages_1q","two_clone_res_collages_2p","two_clone_res_collages_16q",
           "low_hypoxia_clone_tree_files","low_hypoxia_clone_trees_segments_files")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "targets/branches (two-clone + low-hypoxia trees)\n")
stopifnot(length(to_inv) > 0)
' || { echo "FATAL: invalidation failed -- aborting"; exit 1; }

Rscript run_targets_bg.R \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  low_hypoxia_clone_tree_files low_hypoxia_clone_trees_segments_files \
  hypoxia_rebuilt_gdrive

echo "=== output tallies ==="
echo "two-clone collages: $(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -Ec '__(1q|2p|16q)_res')"
echo "low-hypoxia clone trees:   $(ls results/clone_trees/*_low_hypoxia_clone_tree.pdf 2>/dev/null | wc -l)"
echo "low-hypoxia segment trees: $(ls results/clone_trees/*_low_hypoxia_segment_tree.pdf 2>/dev/null | wc -l)"

# Slack ping on completion (best-effort).
if [ -f "$HOME/.slack_webhook" ]; then
  N=$(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -Ec '__(1q|2p|16q)_res')
  T=$(ls results/clone_trees/*_low_hypoxia_*_tree.pdf 2>/dev/null | wc -l)
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"two_clone+trees ($SLURM_JOB_ID) done: $N two-clone collages (SCNA-labeled bars), $T low-hypoxia trees refreshed\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
