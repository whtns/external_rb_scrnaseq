#!/usr/bin/env Rscript
# Pre-register source tarballs in pkgcache by downloading one at a time.
# Run on login node (has internet). After this, sbatch can run on any compute node.
# Must match the sbatch environment: source packages only, x86_64-pc-linux-gnu platform.

Sys.setenv(PKGCACHE_PLATFORMS = "x86_64-pc-linux-gnu")

lockfile <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_source.lock"
pkgs <- jsonlite::fromJSON(lockfile, simplifyVector = TRUE)

# CRAN/Bioc packages by name, GitHub packages by full ref
github_refs <- c(
  "github::mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c",
  "github::cobriniklab/seuratTools@02d63235ca641a491acc756817c432aa686fc73e",
  "github::CCICB/TCGAgistic@adf0782554c2d7916169cf81f69d19247047b7eb",
  "github::mahmoudibrahim/genesorteR@08e234c06bb3bb70ecd8a7e69cecca3cbaa92428",
  "github::immunogenomics/presto@a24772a135c7895a8183b007376050556c60a05b"
)
github_pkgs <- c("SeuratDisk", "seuratTools", "TCGAgistic", "genesorteR", "presto")

cran_names <- pkgs$packages$package[!pkgs$packages$package %in% github_pkgs]
all_refs <- c(cran_names, github_refs)

cat("Downloading", length(all_refs), "packages (source, x86_64-pc-linux-gnu)...\n")

failed <- character(0)
for (i in seq_along(all_refs)) {
  ref <- all_refs[i]
  tryCatch({
    suppressMessages(pak:::pkg_download(ref, dest_dir = tempdir()))
    cat(sprintf("[%d/%d] OK: %s\n", i, length(all_refs), ref))
  }, error = function(e) {
    cat(sprintf("[%d/%d] SKIP: %s\n", i, length(all_refs), ref))
    failed <<- c(failed, ref)
  })
}

cat("\nDone.", length(all_refs) - length(failed), "succeeded,", length(failed), "failed.\n")
if (length(failed)) cat("Failed:", paste(failed, collapse = ", "), "\n")
