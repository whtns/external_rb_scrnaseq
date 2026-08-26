#!/bin/bash
#SBATCH --job-name=collages_recompute
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/collages_recompute_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/collages_recompute_%j.err

# Rebuild BOTH collage sets with the PCA-recompute mitigation, then rclone both:
#  1. heatmap_collages_low_hypoxia_by_res -- low-hypoxia sweep, CAPPED at 0.2..0.8,
#     recompute_pca=TRUE (fresh PCA on the low object, re-cluster each resolution).
#  2. two_clone_res_collages_{1q,2p,16q} -- two-clone SCNA collages, 0.2..1.2,
#     recompute_pca=TRUE: PCA recomputed + clustered on the FULL low object BEFORE
#     the two-clone subset; NO clustering after the subset.
# Then hypoxia_rebuilt_gdrive syncs everything to Drive, grouping the outputs
# under low_hypoxia_by_res/ and two_clone_scna_collage/ (paper_retained/other).

module load r/4.4.1
module load poppler/23.04.0
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers. NAMESPACE is HAND-MAINTAINED; verify the exported
# surface AND that both collage builders carry the recompute_pca arg (guards
# against a stale install silently ignoring the mitigation).
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
ex <- getNamespaceExports("numbatHelpers")
need <- c("plot_hypoxia_low_res_collages","plot_scna_two_clone_res_collages","plot_seu_marker_heatmap")
miss <- setdiff(need, ex); if (length(miss)) stop("NOT EXPORTED: ", paste(miss, collapse=", "))
for (f in c("plot_hypoxia_low_res_collages","plot_scna_two_clone_res_collages"))
  if (!"recompute_pca" %in% names(formals(getExportedValue("numbatHelpers", f))))
    stop(f, " lacks recompute_pca arg -- stale install")
cat("numbatHelpers installed; recompute_pca present on both builders\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed -- aborting"; exit 1; }

# Both target commands changed (capped/recompute_pca) -> outdated, but invalidate
# base targets + dynamic branches explicitly so every branch re-runs.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("heatmap_collages_low_hypoxia_by_res",
           "two_clone_res_collages_1q","two_clone_res_collages_2p","two_clone_res_collages_16q")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "collage targets/branches\n")
' || { echo "FATAL: invalidation failed -- aborting"; exit 1; }

Rscript run_targets_bg.R \
  heatmap_collages_low_hypoxia_by_res \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  hypoxia_rebuilt_gdrive

# Slack ping on completion (best-effort).
if [ -f "$HOME/.slack_webhook" ]; then
  LOW=$(ls results/*_hypoxia_low_seu.rds__filtered_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | wc -l)
  TWO=$(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -E "__(1q|2p|16q)_res" | wc -l)
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"collages_recompute ($SLURM_JOB_ID) done: $LOW low-hypoxia + $TWO two-clone collages on disk, synced to Drive\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
