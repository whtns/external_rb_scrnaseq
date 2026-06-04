#!/usr/bin/env Rscript
# Retry the 10 packages that failed due to libRlapack.so/libRblas.so
# Run with: export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH

options(
  repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/focal/latest",
    BioCsoft = "https://bioconductor.org/packages/3.20/bioc",
    BioCann = "https://bioconductor.org/packages/3.20/data/annotation",
    BioCexp = "https://bioconductor.org/packages/3.20/data/experiment"
  ),
  pak.num_workers = 2,
  Ncpus = 2
)

lib <- .libPaths()[1]
cat("Installing to:", lib, "\n")
cat("LD_LIBRARY_PATH:", Sys.getenv("LD_LIBRARY_PATH"), "\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

# Verify Matrix loads (the key test)
tryCatch({
  library(Matrix)
  cat("Matrix OK:", as.character(packageVersion("Matrix")), "\n\n")
}, error = function(e) {
  stop("Matrix still broken: ", conditionMessage(e))
})

refs <- c(
  "ggtangle",
  "forecast",
  "DESeq2",
  "speckle",
  "plyranges",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "maftools",
  "plotgardener",
  "mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c",
  "CCICB/TCGAgistic@adf0782554c2d7916169cf81f69d19247047b7eb"
)

cat(sprintf("Retrying %d packages\n\n", length(refs)))

for (i in seq_along(refs)) {
  ref <- refs[[i]]
  pkg <- sub("@.*", "", sub(".*/", "", ref))
  cat(sprintf("[%d/%d] Installing %s ... ", i, length(refs), pkg))
  tryCatch({
    pak::pak(ref, lib = lib, ask = FALSE, upgrade = FALSE)
    cat("OK\n")
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
  })
  gc()
}

cat("\nFinal check:\n")
check_pkgs <- c("ggtangle", "forecast", "DESeq2", "speckle", "plyranges",
                "maftools", "plotgardener", "SeuratDisk", "TCGAgistic")
installed <- rownames(installed.packages(lib.loc = lib))
for (p in check_pkgs) {
  cat(sprintf("  %s: %s\n", p, if (p %in% installed) "installed" else "MISSING"))
}
