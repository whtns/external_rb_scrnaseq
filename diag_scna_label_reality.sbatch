#!/bin/bash
#SBATCH --job-name=scna_label_check
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/scna_label_check_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/scna_label_check_%j.err

# Read-only. Does NOT touch _targets_r431, so it is safe to run alongside a
# tar_make(). Settles one question the audit cannot answer from TSVs alone:
# what scna values do the saved Seurat objects ACTUALLY carry?
#
# filter_cluster_save_seu() does large_clone_simplifications[[sample_id]], but
# the target it is fed (large_clone_simplifications_per_sample) returns a bare
# label->seg list, not a one-element list keyed by sample_id the way its
# siblings (large_filter_expressions_per_sample) do. If that indexing misses,
# the key is empty for every sample and EVERY cell is labelled "" = diploid.

set -euo pipefail
module purge
module load gcc/13.3.0   # r/4.4.1 is nested under gcc in the Lmod hierarchy
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
suppressPackageStartupMessages(library(Seurat))
for (sid in c("SRX22868103", "SRX10264525", "SRX10264519")) {
  f <- sprintf("output/seurat/%s_hypoxia_low_seu.rds", sid)
  if (!file.exists(f)) { cat(sid, "missing\n"); next }
  seu <- readRDS(f)
  md <- seu@meta.data
  cat("\n=====", sid, " cells:", ncol(seu), "=====\n")
  if (!"scna" %in% colnames(md)) {
    cat("  NO scna COLUMN -- assemble_diploid_seu() would see NULL here\n")
  } else {
    v <- md$scna; v[is.na(v)] <- "<NA>"; v[v == ""] <- "<empty>"
    print(sort(table(v), decreasing = TRUE))
    cat("  cells passing the diploid mask:",
        sum(is.na(md$scna) | md$scna == ""), "of", nrow(md), "\n")
  }
  if ("GT_opt" %in% colnames(md)) {
    cat("  GT_opt x scna cross-tab:\n")
    print(table(ifelse(is.na(md$GT_opt) | md$GT_opt == "", "GT empty", "GT has SCNA"),
                ifelse(is.na(md$scna)  | md$scna  == "", "scna diploid", "scna labelled")))
  }
  rm(seu, md); invisible(gc(FALSE))
}
'
echo "### exit $?"
