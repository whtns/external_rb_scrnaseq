#!/bin/bash
#SBATCH --job-name=diffex_and_two_clone
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=96G
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diffex_and_two_clone_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/diffex_and_two_clone_%j.err

# ONE job, sequential passes: only one tar_make may touch _targets_r431 at a time.
#
# A) Two-clone SCNA collages now CLUSTER ON THE INCLUDED CELLS.
#    plot_scna_two_clone_res_collages() subsets to the two compared clones first,
#    then recomputes PCA + SNN + FindClusters on that subset, and no longer passes
#    display_cells. This REPLACES the previous lock-to-full zoom, so cluster ids
#    are no longer comparable with the full-population collages or across SCNAs.
#
# B) Three diffex fixes, all verified by src/smoke_diffex_error_fixes.R:
#    1. find_diffex_clones() derived its config key with a "_filtered_seu.*"-only
#       pattern, which does nothing to a *_hypoxia_low_seu.rds path -- so
#       large_clone_comparisons[[sample_id]] was NULL and the function returned an
#       EMPTY list for all 33 samples with no error. That is why every straight
#       (non-per-cluster) by_chr workbook was an empty placeholder. Now uses
#       .seu_sample_id(), which strips the object-type suffix and keeps _branch_N.
#    2. make_volcano_diffex_clones() deduped on (clone_comparison, cluster, chr,
#       symbol) but volcano_plot_clone_clusters() then calls
#       column_to_rownames("symbol"). Genes annotated to a patch contig as well as
#       their real chromosome (ATMIN -> chr16 AND HG405_PATCH) survived as two rows
#       -> "duplicate 'row.names' are not allowed". Dedupe key now matches the
#       rowname key, canonical chromosome retained.
#    3. tally_kooi_candidates() died with "object 'kooi_region' not found" on the
#       empty placeholder workbook (write_xlsx of an empty list -> a bare Sheet1
#       with no columns). Now skips column-less sheets and returns empty.

module load r/4.4.1
module load poppler/23.04.0
module load rclone
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install for the crew workers, then assert the fixes are actually in the
# INSTALLED package -- a stale install silently reproduces every bug above.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
ns <- asNamespace("numbatHelpers")
ex <- getNamespaceExports("numbatHelpers")
need <- c("plot_scna_two_clone_res_collages","plot_seu_marker_heatmap",
          "find_diffex_clones","tally_kooi_candidates","make_volcano_diffex_clones")
miss <- setdiff(need, ex); if (length(miss)) stop("NOT EXPORTED: ", paste(miss, collapse=", "))
for (f in c(".seu_sample_id",".canonical_chr_first",".read_kooi_candidates",".read_kooi_sheets"))
  if (!exists(f, ns)) stop("STALE install: ", f, " absent")
# functional assertions, not just presence
sid <- get(".seu_sample_id", ns)
stopifnot(sid("output/seurat/SRX10031194_hypoxia_low_seu.rds") == "SRX10031194")
stopifnot(sid("output/seurat/SRX10264523_branch_5_filtered_seu.rds") == "SRX10264523_branch_5")
k <- get(".read_kooi_candidates", ns)()
stopifnot(!any(duplicated(k$symbol[k$symbol == "TDP2"])))
cat("numbatHelpers installed; sample-id, kooi, and chr fixes verified\n")' \
  || { echo "FATAL: install/verify failed -- aborting before tar_make"; exit 1; }

# targets does not hash package-internal functions, and these targets map over
# path-STRING targets (seus_low_hypoxia), so nothing here auto-invalidates.
# Invalidate explicitly:
#   - two_clone_res_collages_*    : re-cluster-on-included-cells behaviour change
#   - {all,cis,trans}_diffex_clones : were cached EMPTY by the sample-id bug
#   - volcano_*, table_*, found_kooi_candidates : consumers of the above
# The *_for_each_cluster families are deliberately NOT invalidated: they key on
# tumor_id (never hit the sample-id bug), already built 33/33 with 0 errors, and
# rebuilding them costs hours for no change.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("two_clone_res_collages_1q","two_clone_res_collages_2p","two_clone_res_collages_16q",
           "all_diffex_clones","cis_diffex_clones","trans_diffex_clones",
           "volcano_all","volcano_cis","volcano_trans",
           "table_all_diffex_clones","table_cis_diffex_clones","table_trans_diffex_clones",
           "found_kooi_candidates")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "targets/branches\n")
stopifnot(length(to_inv) > 0)
' || { echo "FATAL: invalidation failed -- aborting"; exit 1; }

echo "=== pass 1: two-clone collages (re-clustered on included cells) + drive sync ==="
Rscript run_targets_bg.R \
  two_clone_res_collages_1q two_clone_res_collages_2p two_clone_res_collages_16q \
  hypoxia_rebuilt_gdrive \
  || echo "WARN: two-clone collage pass failed"

# found_kooi_candidates reads the by_chr xlsx as HARDCODED string paths, not target
# dependencies, so in a single tar_make it is dispatched immediately -- before the
# tables have written those files -- and under error="stop" its failure aborts the
# whole run. Build the tables/volcanos first, then the kooi tally.
echo "=== pass 2: diffex tables + volcanos ==="
Rscript run_targets_bg.R \
  volcano_all volcano_cis volcano_trans \
  table_all_diffex_clones table_cis_diffex_clones table_trans_diffex_clones \
  || echo "WARN: diffex tables/volcanos pass failed"

echo "=== pass 3: found_kooi_candidates ==="
Rscript run_targets_bg.R found_kooi_candidates || echo "WARN: found_kooi_candidates failed"

echo "=== two-clone collage tally ==="
echo "two-clone collages: $(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -Ec '__(1q|2p|16q)_res')"

echo "=== diffex deliverables ==="
for f in results/diffex_bw_clones_per_cluster_all.xlsx \
         results/diffex_bw_clones_per_cluster_large_in_segment.xlsx \
         results/diffex_bw_clones_per_cluster_large_out_of_segment.xlsx \
         results/diffex_bw_clones_per_cluster_all.pdf \
         results/diffex_bw_clones_per_cluster_large_in_segment.pdf \
         results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf \
         results/diffex_bw_clones_all.xlsx \
         results/diffex_bw_clones_all_by_chr.xlsx \
         results/diffex_bw_clones_large_in_segment_by_chr.xlsx \
         results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx; do
  if [ -f "$f" ]; then echo "  OK   $(stat -c '%y' "$f" | cut -d. -f1)  $(du -h "$f" | cut -f1)  $f"
  else echo "  MISSING  $f"; fi
done

# The straight by_chr workbooks were 4937-byte empty placeholders before the
# sample-id fix; report whether they now carry real sheets.
Rscript -e '
suppressPackageStartupMessages(library(readxl))
for (f in c("results/diffex_bw_clones_all_by_chr.xlsx",
            "results/diffex_bw_clones_large_in_segment_by_chr.xlsx",
            "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx")) {
  if (!file.exists(f)) { cat("  MISSING", f, "\n"); next }
  s <- excel_sheets(f)
  cat(sprintf("  %-62s sheets: %s\n", basename(f), paste(head(s, 12), collapse=",")))
}' || true

echo "=== errored targets ==="
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
m <- tar_meta(complete_only=FALSE)
for (f in c("all_diffex_clones","cis_diffex_clones","trans_diffex_clones",
            "all_diffex_clones_for_each_cluster","cis_diffex_clones_for_each_cluster",
            "trans_diffex_clones_for_each_cluster",
            "two_clone_res_collages_1q","two_clone_res_collages_2p","two_clone_res_collages_16q")) {
  r <- m[grepl(paste0("^",f,"(_[0-9a-f]+)?$"), m$name), ]
  cat(sprintf("  %-40s branches: %2d | errored: %d\n", f, nrow(r), sum(!is.na(r$error))))
}
e <- m[!is.na(m$error) & grepl("^(volcano_|found_kooi|table_)", m$name), c("name","error")]
if (nrow(e)) for (i in seq_len(nrow(e))) cat("  ERR", e$name[i], ":", substr(e$error[i],1,160), "\n") else
  cat("  no volcano/table/kooi errors\n")
' || true

if [ -f "$HOME/.slack_webhook" ]; then
  N=$(ls results/*_hypoxia_low_seu.rds__*_res*_heatmap_phase_scatter_patchwork.pdf 2>/dev/null | grep -Ec '__(1q|2p|16q)_res')
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"diffex+two_clone ($SLURM_JOB_ID) done: $N two-clone collages (re-clustered on included cells), diffex rebuilt\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
