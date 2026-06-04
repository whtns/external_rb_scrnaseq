#!/usr/bin/env Rscript
# Install wrong-version packages directly via install.packages() — no pak, no cache issues.
# Run on a node WITH internet (login node or OnDemand).
# Usage: Rscript install_missing.R [--dry-run]

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.4")
lockfile <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_source.lock"
logfile <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj/install_missing.log"

cat("Library:", lib, "\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

pkgs <- jsonlite::fromJSON(lockfile, simplifyVector = TRUE)[["packages"]]
github_types <- c("github")

# Packages to skip:
# - Matrix/KernSmooth: base R packages at newer versions (can't downgrade)
# - ggtree: not used in pipeline; 3.14.0 incompatible with ggplot2 4.0 (check_linewidth removed)
skip_packages <- c("Matrix", "KernSmooth", "ggtree")

# Set env vars needed for packages that use spack-provided libraries.
# CURL_ROOT/BZIP2_ROOT: Rsamtools/rtracklayer src/Makevars expand -L$(X_ROOT)/lib
# CAIRO_CFLAGS: cairo.h is at include/cairo/cairo.h, so include/cairo/ needed
# PKG_CONFIG_PATH: already set by module load cairo/fontconfig/freetype (passed in via env)
spack_gcc  <- "/apps/spack/2406/apps/linux-rocky8-x86_64_v3/gcc-13.3.0"
curl_root  <- Sys.getenv("CURL_ROOT",  paste0(spack_gcc, "/curl-8.8.0-sqn4mlx"))
bzip2_root <- Sys.getenv("BZIP2_ROOT", paste0(spack_gcc, "/bzip2-1.0.8-lpw5x7j"))
cairo_root <- Sys.getenv("CAIRO_ROOT", paste0(spack_gcc, "/cairo-1.16.0-iva3nxi"))
Sys.setenv(
  CURL_ROOT  = curl_root,
  BZIP2_ROOT = bzip2_root,
  CAIRO_CFLAGS = paste0("-I", cairo_root, "/include/cairo"),
  CAIRO_LIBS   = paste0("-L", cairo_root, "/lib -lcairo")
)

hdf5_h5cc <- "/apps/spack/2406/apps/linux-rocky8-x86_64_v3/gcc-13.3.0/hdf5-1.12.3-s4big4c/bin/h5cc"

# Build list of wrong-version packages (version-normalized comparison)
get_needs_install <- function() {
  installed <- installed.packages(lib.loc = lib)[, c("Package", "Version")]
  needs <- pkgs[mapply(function(nm, ver) {
    if (nm %in% skip_packages) return(FALSE)
    idx <- which(installed[, "Package"] == nm)
    if (length(idx) == 0) return(TRUE)
    inst_ver <- installed[idx, "Version"]
    !tryCatch(
      numeric_version(inst_ver) == numeric_version(gsub("-", ".", ver)),
      error = function(e) inst_ver == ver
    )
  }, pkgs$package, pkgs$version), ]
  needs
}

needs <- get_needs_install()
cat("Packages needing install/update:", nrow(needs), "\n\n")

if (dry_run) {
  cat("DRY RUN — packages that would be installed:\n")
  print(needs[, c("package", "version", "type", "sources")])
  quit(status = 0)
}

ver_match <- function(inst_ver, wanted_ver) {
  # Normalize both sides: package_version converts "-" to "." (e.g. 1.98-1.18 == 1.98.1.18)
  tryCatch(
    numeric_version(inst_ver) == numeric_version(gsub("-", ".", wanted_ver)),
    error = function(e) inst_ver == wanted_ver
  )
}

install_one_cran <- function(nm, ver, url, idx, total) {
  if (is.list(url)) url <- url[[1]]
  cat(sprintf("[%d/%d] Installing %s %s\n  from %s\n", idx, total, nm, ver, url))
  # Package-specific configure arguments
  conf_args <- switch(nm,
    hdf5r = paste0("--with-hdf5=", hdf5_h5cc),
    NULL
  )
  tryCatch({
    if (!is.null(conf_args)) {
      install.packages(url, repos=NULL, type="source", lib=lib, dependencies=FALSE,
                       quiet=FALSE, configure.args=conf_args)
    } else {
      install.packages(url, repos=NULL, type="source", lib=lib, dependencies=FALSE, quiet=FALSE)
    }
    inst <- tryCatch(as.character(packageVersion(nm, lib.loc=lib)), error=function(e) NULL)
    if (is.null(inst) || !ver_match(inst, ver))
      stop(sprintf("version mismatch after install: got %s, wanted %s",
                   if (is.null(inst)) "NA" else inst, ver))
    cat(sprintf("  ✔ %s %s installed\n", nm, ver))
    TRUE
  }, error=function(e) {
    msg <- conditionMessage(e)
    # Bioc packages not on CRAN — fall back to BiocManager (only on download/404 errors, not version mismatches)
    if (grepl("cannot open|404|does not exist|no packages available|download.file", msg, ignore.case=TRUE)) {
      cat(sprintf("  ↺ Trying BiocManager for %s %s\n", nm, ver))
      bioc_ok <- tryCatch({
        BiocManager::install(nm, version=BiocManager::version(), lib=lib,
                             update=FALSE, ask=FALSE, quiet=FALSE)
        inst2 <- tryCatch(as.character(packageVersion(nm, lib.loc=lib)), error=function(e) NULL)
        if (!is.null(inst2) && ver_match(inst2, ver)) {
          cat(sprintf("  ✔ %s %s installed via BiocManager\n", nm, ver))
          TRUE
        } else {
          cat(sprintf("  ✖ BiocManager got %s, wanted %s\n", inst2, ver))
          FALSE
        }
      }, error=function(e2) {
        cat(sprintf("  ✖ BiocManager FAILED: %s\n", conditionMessage(e2)))
        FALSE
      })
      return(bioc_ok)
    }
    cat(sprintf("  ✖ FAILED: %s\n", msg))
    FALSE
  })
}

github_refs <- list(
  SeuratDisk      = "mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c",
  seuratTools     = "cobriniklab/seuratTools@02d63235ca641a491acc756817c432aa686fc73e",
  TCGAgistic      = "CCICB/TCGAgistic@adf0782554c2d7916169cf81f69d19247047b7eb",
  genesorteR      = "mahmoudibrahim/genesorteR@08e234c06bb3bb70ecd8a7e69cecca3cbaa92428",
  presto          = "immunogenomics/presto@a24772a135c7895a8183b007376050556c60a05b",
  numbat          = "whtns/numbat@b9354cd425008d1ac228abd8266b3ae05f07e8eb",
  ComplexHeatmap  = "jokergoo/ComplexHeatmap@ad11b26afe8c63dc05b20035503352ca740f602e",
  DataEditR       = "DillonHammill/DataEditR@76985bd2c39ecac5d6f56dbb1f8ada85bcd5d85f"
)

# ── Multi-pass loop: retry until no progress ──────────────────────────────────
max_passes <- 5
for (pass in seq_len(max_passes)) {
  needs <- get_needs_install()
  if (nrow(needs) == 0) break
  cat(sprintf("\n══ Pass %d: %d packages to install ══\n\n", pass, nrow(needs)))

  github_pkgs <- needs[needs$type %in% github_types, ]
  cran_pkgs   <- needs[!needs$type %in% github_types, ]
  cat("CRAN/Bioc:", nrow(cran_pkgs), "  GitHub:", nrow(github_pkgs), "\n\n")

  installed_this_pass <- 0

  for (i in seq_len(nrow(cran_pkgs))) {
    ok <- install_one_cran(cran_pkgs$package[i], cran_pkgs$version[i],
                           cran_pkgs$sources[i], i, nrow(cran_pkgs))
    if (ok) installed_this_pass <- installed_this_pass + 1
  }

  if (nrow(github_pkgs) > 0) {
    cat("\n── GitHub packages ──\n")
    for (nm in github_pkgs$package) {
      ref <- github_refs[[nm]]
      if (is.null(ref)) { cat(sprintf("  ? %s: no GitHub ref, skipping\n", nm)); next }
      cat(sprintf("  Installing %s from %s\n", nm, ref))
      tryCatch({
        remotes::install_github(ref, lib=lib, upgrade="never", quiet=FALSE)
        cat(sprintf("  ✔ %s installed\n", nm))
        installed_this_pass <- installed_this_pass + 1
      }, error=function(e) cat(sprintf("  ✖ FAILED: %s\n", conditionMessage(e))))
    }
  }

  cat(sprintf("\nPass %d done: %d/%d installed\n", pass, installed_this_pass, nrow(needs)))
  if (installed_this_pass == 0) {
    cat("No progress in this pass — stopping to avoid infinite loop\n")
    break
  }
}

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n══ Final Summary ══\n")
remaining <- get_needs_install()
cat("Still wrong-version after all passes:", nrow(remaining), "\n")
if (nrow(remaining) > 0) {
  cat("Remaining:\n")
  print(remaining[, c("package", "version")])
}
if (nrow(remaining) == 0) cat("All packages at correct versions!\n")
cat("Done:", format(Sys.time()), "\n")
