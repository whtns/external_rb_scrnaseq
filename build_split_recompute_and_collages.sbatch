#!/bin/bash
#SBATCH --job-name=split_recompute_collages
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/split_recompute_collages_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/split_recompute_collages_%j.err

# Rebuild the hypoxia split so each *_hypoxia_low_seu.rds is created with a PCA
# recomputed on the low-hypoxia cells (subset_seu_by_expression recompute_pca=TRUE),
# then rebuild BOTH collage sets (which now read the fresh persisted clusters) and
# rclone everything to Drive.
#
# hypoxia_partition_paths carries cue=tar_cue(command=FALSE) (deliberately frozen),
# so the subset_seu_by_expression change does NOT auto-invalidate it -- invalidate
# it explicitly. That reruns the split for ALL samples (regenerating low+high
# objects; the low/high MEMBERSHIP is unchanged -- only the low object's embedding
# + clusters change) and cascades to seus_low_hypoxia, hypoxia_seus_*, and the
# collages.

module load r/4.4.1
module load poppler/23.04.0
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers; verify the recompute_pca arg reached the installed
# package (NAMESPACE is hand-maintained; guard against a stale install).
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
ns <- asNamespace("numbatHelpers")
if (!"recompute_pca" %in% names(formals(get("subset_seu_by_expression", ns))))
  stop("subset_seu_by_expression lacks recompute_pca -- stale install")
ex <- getNamespaceExports("numbatHelpers")
need <- c("plot_hypoxia_low_res_collages","plot_scna_two_clone_res_collages")
miss <- setdiff(need, ex); if (length(miss)) stop("NOT EXPORTED: ", paste(miss, collapse=", "))
cat("numbatHelpers installed; subset_seu_by_expression has recompute_pca; collage builders exported\n")' \
  || { echo "FATAL: install/verify failed -- aborting"; exit 1; }

# Invalidate the frozen split (base + branches) so it reruns with recompute_pca,
# plus the collage targets so every branch re-renders from the fresh low objects.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("hypoxia_partition_paths",
           "heatmap_collages_low_hypoxia_by_res",
           "two_clone_res_collages_1q","two_clone_res_collages_2p","two_clone_res_collages_16q")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "targets/branches (split + collages)\n")
stopifnot(length(to_inv) > 0)
' || { echo "FATAL: invalidation failed -- aborting"; exit 1; }

# Build: split (regenerates low objects) -> collages -> Drive sync. targets pulls
# hypoxia_partition_paths + seus_low_hypoxia + hypoxia_seus_* as upstream deps.
Rscript run_targets_bg.R \
  heatmap_collages_low_hypoxia_by_res \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  hypoxia_rebuilt_gdrive

# Report any low objects that came back NULL (error="null" caches silent NULLs).
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
low <- tryCatch(unlist(tar_read(seus_low_hypoxia)), error=function(e) character(0))
n_na <- sum(is.na(low) | low=="")
cat("seus_low_hypoxia:", length(low), "paths;", n_na, "NA/empty (failed split branches)\n")
' || true

if [ -f "$HOME/.slack_webhook" ]; then
  LOW=$(ls results/*_hypoxia_low_seu.rds__filtered_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | wc -l)
  TWO=$(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -E "__(1q|2p|16q)_res" | wc -l)
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"split_recompute_collages ($SLURM_JOB_ID) done: $LOW low-hypoxia + $TWO two-clone collages on disk, synced to Drive\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
