#!/usr/bin/env Rscript
# Install missing packages from pak.lock one at a time to avoid OOM

options(
  repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/focal/latest",
    BioCsoft = "https://bioconductor.org/packages/3.20/bioc",
    BioCann = "https://bioconductor.org/packages/3.20/data/annotation",
    BioCexp = "https://bioconductor.org/packages/3.20/data/experiment"
  )
)

# Limit pak workers to avoid OOM
options(
  pak.num_workers = 2,
  Ncpus = 2
)

lib <- .libPaths()[1]
cat("Installing to:", lib, "\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

lockfile <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak.lock"

lock_data <- jsonlite::fromJSON(lockfile, simplifyVector = FALSE)

installed <- rownames(installed.packages(lib.loc = lib))

missing_pkgs <- Filter(function(p) {
  p$type != "installed" && !(p$package %in% installed)
}, lock_data$packages)

cat(sprintf("Missing packages: %d\n\n", length(missing_pkgs)))

refs <- sapply(missing_pkgs, function(p) p$ref)
cat("Refs to install:\n")
cat(paste(refs, collapse = "\n"), "\n\n")

# Install one at a time with gc() between each
for (i in seq_along(refs)) {
  ref <- refs[[i]]
  pkg <- missing_pkgs[[i]]$package
  cat(sprintf("[%d/%d] Installing %s (ref: %s) ... ", i, length(refs), pkg, ref))
  tryCatch({
    pak::pak(ref, lib = lib, ask = FALSE, upgrade = FALSE)
    cat("OK\n")
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
  })
  gc()
  Sys.sleep(1)
}

cat("\nDone. Re-checking installed:\n")
installed_after <- rownames(installed.packages(lib.loc = lib))
still_missing <- sapply(missing_pkgs, function(p) p$package)
still_missing <- still_missing[!(still_missing %in% installed_after)]
if (length(still_missing) == 0) {
  cat("All packages successfully installed!\n")
} else {
  cat("Still missing:", paste(still_missing, collapse = ", "), "\n")
}
