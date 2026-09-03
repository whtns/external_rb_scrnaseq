#!/bin/bash
#SBATCH --job-name=validate_dag
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/validate_dag_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/validate_dag_%j.err

# Validate the pipeline definition after the #41/#42/#43/#44 edits.
#
# tar_manifest() parses _targets.R and every R/pipeline_targets_*.R and resolves
# the dependency graph WITHOUT building anything. It is the cheapest way to
# catch a dangling reference to a target I retired (filtered_numbat_pdfs,
# low_hypoxia_numbat_pdfs, ks_diploid_seu) or a typo in a new one.
#
# Nothing is built and nothing is invalidated, so this does NOT count as a
# tar_make() against _targets_r431 and is safe to run on its own.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# install, not load_all: the pipeline definition calls library(numbatHelpers)
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers",
                                         force = TRUE, upgrade = "never", quiet = TRUE))
cat("numbatHelpers installed\n")' || { echo "FATAL: install failed"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store = "_targets_r431")

m <- tar_manifest()
cat("targets defined:", nrow(m), "\n\n")

# the ones this work added or changed
new <- c("numbat_bulk_clone_pdfs", "unfiltered_numbat_bulk_clones",
         "filtered_numbat_bulk_clones", "low_hypoxia_numbat_bulk_clones",
         "unfiltered_numbat_expression", "filtered_numbat_expression",
         "low_hypoxia_numbat_expression",
         "diploid_panel_samples_file", "diploid_panel_sample_ids",
         "diploid_seu", "diploid_seu_composition")
cat("=== new / changed targets present ===\n")
for (n in new) cat(sprintf("  %-32s %s\n", n, n %in% m$name))

# the ones this work retired
gone <- c("filtered_numbat_pdfs", "low_hypoxia_numbat_pdfs", "ks_diploid_seu")
cat("\n=== retired targets absent ===\n")
for (n in gone) cat(sprintf("  %-32s %s\n", n, !(n %in% m$name)))

ok <- all(new %in% m$name) && !any(gone %in% m$name)

# A dangling dependency would make the network fail to resolve.
net <- tar_network(targets_only = TRUE)
missing <- setdiff(unique(c(net$edges$from, net$edges$to)), m$name)
cat("\n=== dangling dependencies ===\n")
if (length(missing)) {
  cat("  ", paste(missing, collapse=", "), "\n"); ok <- FALSE
} else {
  cat("  none\n")
}

cat("\n", if (ok) "DAG OK\n" else "DAG PROBLEM\n", sep="")
if (!ok) quit(status = 1)
'
