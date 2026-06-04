#!/usr/bin/env Rscript
# Fix remaining 2 packages:
#   forecast: needs urca from source (Rocky 8 binary requires GLIBC_2.29)
#   SeuratDisk: needs libhdf5.so.200 in LD_LIBRARY_PATH (now in ~/lib)

lib <- .libPaths()[1]
cat("LD_LIBRARY_PATH:", Sys.getenv("LD_LIBRARY_PATH"), "\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

# 1. Reinstall urca from source (force compile on this system's glibc 2.28)
cat("[1] Reinstalling urca from source...\n")
tryCatch({
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/urca/urca_1.3-3.tar.gz",
    repos = NULL, type = "source", lib = lib
  )
  # Try urca 1.3-4 (match pak.lock version)
  install.packages("urca", repos = "https://cloud.r-project.org",
                   type = "source", lib = lib)
  cat("urca OK\n\n")
}, error = function(e) cat("urca FAILED:", conditionMessage(e), "\n\n"))

# 2. Now install forecast (depends on urca)
cat("[2] Installing forecast...\n")
tryCatch({
  pak::pak("forecast", lib = lib, ask = FALSE, upgrade = FALSE)
  cat("forecast OK\n\n")
}, error = function(e) {
  cat("forecast pak FAILED, trying install.packages from source...\n")
  tryCatch({
    install.packages("forecast", repos = "https://cloud.r-project.org",
                     type = "source", lib = lib)
    cat("forecast OK (source)\n\n")
  }, error = function(e2) cat("forecast FAILED:", conditionMessage(e2), "\n\n"))
})

# 3. SeuratDisk (needs libhdf5.so.200)
cat("[3] Installing SeuratDisk...\n")
tryCatch({
  pak::pak("mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c",
           lib = lib, ask = FALSE, upgrade = FALSE)
  cat("SeuratDisk OK\n\n")
}, error = function(e) {
  cat("SeuratDisk FAILED:", conditionMessage(e), "\n\n")
})

# Final summary
cat("=== Final check ===\n")
installed <- rownames(installed.packages(lib.loc = lib))
for (p in c("forecast", "urca", "SeuratDisk")) {
  status <- if (p %in% installed) "installed" else "MISSING"
  cat(sprintf("  %s: %s\n", p, status))
}
