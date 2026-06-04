#!/usr/bin/env Rscript
# Full pak restore from lockfile with limited workers to avoid OOM
# Run with: bash --login -c 'module load r/4.4.1 2>/dev/null; Rscript ...'

options(
  pak.num_workers = 1,
  Ncpus = 1,
  # Build stringi with bundled ICU so it doesn't link to libicui18n.so.66,
  # which is not available on Rocky Linux 8 (ubuntu binary links to ICU 66).
  configure.args = c(stringi = "--disable-pkg-config")
)

lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.4")
cat("Library:", lib, "\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

lockfile <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_source.lock"
cat("Running pak::lockfile_install from:", lockfile, "\n\n")

tryCatch({
  pak::lockfile_install(lockfile, lib = lib)
  cat("\nlockfile_install completed successfully\n")
}, error = function(e) {
  cat("\nlockfile_install FAILED:", conditionMessage(e), "\n")
  quit(status = 1)
})

cat("Done:", format(Sys.time()), "\n")
