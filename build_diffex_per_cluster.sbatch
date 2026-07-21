#!/bin/bash
#SBATCH --job-name=diffex_per_cluster
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=96G
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diffex_per_cluster_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diffex_per_cluster_%j.err

# Rebuild the between-clone, per-cluster, per-SCNA differential expression.
#
# find_diffex_bw_clones_for_each_cluster(seus_low_hypoxia, ..., large_clone_comparisons,
# cluster_orders, location=all/cis/out_of_segment) keys clusters on SCT_snn_res.0.6
# and the per-SCNA comparisons come from large_clone_comparisons. The store held two
# STALE errors (both already fixed in committed source, targets just never reran):
#   - cis "`x` must be a vector"  -> guarded in find_diffex_bw_clones_for_each_cluster
#     (derives `clusters` from group.by when the low object lacks it)
#   - table_all_diffex_clones NA-path crash -> .drop_missing_paths() in
#     tabulate_diffex_clones / make_volcano_diffex_clones
# and are stale vs the 2026-07-17 fresh-PCA low-object clusters. These targets map
# over seus_low_hypoxia (path strings), so they never auto-invalidated -> invalidate
# explicitly and rebuild on the CURRENT clusters.

module load r/4.4.1
module load poppler/23.04.0
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for crew workers; verify the two guards are present in the installed pkg.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
ns <- asNamespace("numbatHelpers")
for (f in c("find_diffex_bw_clones_for_each_cluster","tabulate_diffex_clones","make_volcano_diffex_clones"))
  if (!exists(f, ns)) stop("missing fn: ", f)
if (!exists(".drop_missing_paths", ns)) stop("STALE install: .drop_missing_paths guard absent")
# .drop_missing_paths must drop NA/empty/nonexistent
g <- get(".drop_missing_paths", ns)
stopifnot(length(g(list(NA_character_, "", "/no/such/file"))) == 0)
cat("numbatHelpers installed; diffex fns + .drop_missing_paths guard present\n")' \
  || { echo "FATAL: install/verify failed -- aborting"; exit 1; }

# Invalidate ONLY the failed downstream consumers (volcano/table/kooi). The
# _for_each_cluster + _clones branches already built cleanly (all:33, cis:33,
# trans:partial) and are cached; reinstalling the package does NOT invalidate them
# (targets does not hash package fns), so they are reused and the tables merely
# finish any remaining trans branches on demand. The volcano/table reruns pick up
# the newly installed make_volcano_diffex_clones / tabulate_diffex_clones fix.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("volcano_all","volcano_cis","volcano_trans",
           "table_all_diffex_clones","table_cis_diffex_clones","table_trans_diffex_clones",
           "found_kooi_candidates")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "diffex consumer targets\n")
stopifnot(length(to_inv) > 0)
' || { echo "FATAL: invalidation failed -- aborting"; exit 1; }

# Build in TWO ordered passes. found_kooi_candidates reads the cis/trans by_chr
# xlsx as HARDCODED string paths (not target deps), so if requested in the same
# tar_make targets dispatches it immediately -- before the tables write those files
# -- and its error aborts the whole run under error="stop". So build the tables +
# volcanos first (these pull the _for_each_cluster + _clones branches), THEN the
# kooi tally once its input xlsx exist.
Rscript run_targets_bg.R \
  volcano_all volcano_cis volcano_trans \
  table_all_diffex_clones table_cis_diffex_clones table_trans_diffex_clones \
  || { echo "FATAL: diffex tables/volcanos build failed"; }

echo "=== pass 2: found_kooi_candidates (needs the cis/trans by_chr xlsx) ==="
Rscript run_targets_bg.R found_kooi_candidates || echo "WARN: found_kooi_candidates failed"

echo "=== diffex-per-cluster deliverables on disk ==="
for f in results/diffex_bw_clones_per_cluster_all.xlsx \
         results/diffex_bw_clones_per_cluster_large_in_segment.xlsx \
         results/diffex_bw_clones_per_cluster_large_out_of_segment.xlsx \
         results/diffex_bw_clones_per_cluster_all.pdf \
         results/diffex_bw_clones_per_cluster_large_in_segment.pdf \
         results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf; do
  if [ -f "$f" ]; then echo "  OK   $(stat -c '%y' "$f" | cut -d. -f1)  $f"; else echo "  MISSING  $f"; fi
done

# Report any diffex branches that came back errored/NA.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
m <- tar_meta(complete_only=FALSE)
for (f in c("all_diffex_clones_for_each_cluster","cis_diffex_clones_for_each_cluster","trans_diffex_clones_for_each_cluster")) {
  r <- m[grepl(paste0("^",f,"_[0-9a-f]+$"), m$name), ]
  cat(sprintf("  %-38s branches: %d | errored: %d\n", f, nrow(r), sum(!is.na(r$error))))
}
' || true

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"diffex_per_cluster ($SLURM_JOB_ID) rebuild done\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
