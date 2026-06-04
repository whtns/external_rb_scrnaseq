seuratTools is present at /project2/cobrinik_1090/rpkgs/seuratTools

use the ARMOR conda env
apptainer is NOT required for the targets pipeline; run targets directly with Rscript
use module load apptainer
use the apptainer container at /project2/cobrinik_1090/external_rb_scrnaseq_proj/pipeline/containers/numbat-pipeline.sif

## R Environment (Rocky Linux 8 / SLURM)
- Load R: `module load r/4.4.1`
- Always set before running R: `export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"`
  - ~/lib has symlinks: libRlapack.so, libRblas.so → openblas; libhdf5.so.200, libhdf5_hl.so.200 → hdf5-1.12.3
- R library: `~/R/x86_64-pc-linux-gnu-library/4.4`
- Project .Rprofile prepends `.R_libs` to .libPaths(); always pass `lib = path.expand("~/R/...")` explicitly to pak

## R Package Installation
- Lockfile: `/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_source.lock`
  - Modified from pak.lock: all binary=FALSE, platform=source; annotables removed; 429 packages total
  - annotables was OOM-killed during install — use AnnotationDbi + org.Hs.eg.db instead
- Restore script: `/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_restore_full.R`
- Submit install job: `sbatch /project2/cobrinik_1090/external_rb_scrnaseq_proj/install_packages.sbatch`
  - Needs 32 GB RAM (sbatch --mem=32G); pak planning phase OOMs at 16 GB
  - Compute nodes (a01-16) may lack internet — pre-download source tarballs first if needed
  - Pre-download: run a download loop on the login/d05-29 node before sbatch

## Known Gotchas
- **stringi ICU issue**: ubuntu binary links to libicui18n.so.66, which is not on Rocky Linux 8.
  Fix: ensure stringi's `target` in pak_source.lock is `src/contrib/stringi_1.8.7.tar.gz`
  (the standard CRAN source path, NOT the ubuntu-20.04 path). Source tarball must be a real
  .tar.gz (12 MB), not an HTML redirect. Download from CRAN if the cached file is tiny (<1 MB).
  pak_restore_full.R also sets `options(configure.args = c(stringi = "--disable-pkg-config"))`
  as a safety net to force bundled ICU if pkg-config finds the wrong version.
- **pak build cache poisoning**: pak caches ubuntu-derived builds as rocky-8.10 artifacts.
  If a package keeps failing with GLIBC_2.29 or ICU errors: delete
  `~/.cache/R/pkgcache/R/pkgcache/pkg/src/contrib/x86_64-pc-linux-gnu-ubuntu-20.04/` and rerun.
- **Fake source tarballs**: PPM sometimes returns an HTML redirect instead of the tarball.
  Check file size — a real source tarball is at least several hundred KB.
  Re-download from CRAN directly: `wget https://cran.r-project.org/src/contrib/PKG_VER.tar.gz`

## Running the Pipeline
- Script: `Rscript /project2/cobrinik_1090/external_rb_scrnaseq_proj/run_targets_bg.R`
- targets store: `_targets_r431`; default target: `figures_and_tables`
- Set LD_LIBRARY_PATH before running (see R Environment above)
