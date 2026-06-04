always write markdown files to the docs directory
seuratTools is present at /project2/cobrinik_1090/rpkgs/seuratTools

use the ARMOR conda env
apptainer is NOT required for the targets pipeline; run targets directly with Rscript
use module load apptainer
use the apptainer container at /project2/cobrinik_1090/external_rb_scrnaseq_proj/pipeline/containers/numbat-pipeline.sif

## R Environment (Rocky Linux 8 / SLURM)
- Load R: `module load r/4.4.1`
- Always set before running R: `export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"`
  - ~/lib has symlinks: libRlapack.so, libRblas.so → openblas; libhdf5.so.200, libhdf5_hl.so.200 → hdf5-1.12.3
- R library: `~/R/x86_64-pc-linux-gnu-library/4.4` — single location for all packages
- `.R_libs` has been removed; `.Rprofile` no longer prepends it. Do not recreate it.

## R Package Installation
- Lockfile: `/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_source.lock`
  - Modified from pak.lock: all binary=FALSE, platform=source; annotables removed; 429 packages total
  - annotables was OOM-killed during install — use AnnotationDbi + org.Hs.eg.db instead
- Restore script: `/project2/cobrinik_1090/external_rb_scrnaseq_proj/pak_restore_full.R`
- Submit install job: `sbatch /project2/cobrinik_1090/external_rb_scrnaseq_proj/install_packages.sbatch`
  - Needs 32 GB RAM (sbatch --mem=32G); pak planning phase OOMs at 16 GB
  - Compute nodes (a01-16) may lack internet — pre-download source tarballs first if needed
  - Pre-download: run a download loop on the login/d05-29 node before sbatch

## Developing numbatHelpers Package

For fast iteration when editing functions in `/project2/cobrinik_1090/rpkgs/numbat_helpers/`:

**Interactive testing (current R session):**
- Use `devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")` to reload in-memory (instant)
- Workflow: edit R file → run load_all() → test interactively
- This avoids the 30-60s rebuild/install cycle for rapid development

**For sbatch jobs:**
- Add `devtools::load_all()` to the sbatch script to load latest code into the new R process:
```bash
#!/bin/bash
#SBATCH ...
module load r/4.4.1
export LD_LIBRARY_PATH="..."
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
Rscript -e "devtools::load_all('/project2/cobrinik_1090/rpkgs/numbat_helpers')" && \
  Rscript run_targets_bg.R target_name
```
- This loads the latest code without reinstalling the package
- Only use `devtools::install_local()` when you need the package available to all processes (e.g., distributed crew workers)
  - Install to the single library: `devtools::install_local('/project2/cobrinik_1090/rpkgs/numbat_helpers', force = TRUE)` (no `lib =` needed; ~/R/ is the only library path)

**Git commit workflow for numbatHelpers:**
- Automated daily commits: A Stop hook runs at the end of each session to commit uncommitted changes with message `Daily update: YYYY-MM-DD`
- Meaningful development commits: When fixing bugs or adding features, commit with descriptive messages using the format: `type(component): description`
  - Types: `fix` (bug fix), `feat` (new feature), `refactor` (code restructuring), `test` (test changes), `docs` (documentation), `chore` (maintenance)
  - Component: the function or file being changed (e.g., `plot_functions_4`, `NAMESPACE`)
  - Examples: `fix(plot_functions_4): handle empty segments correctly`, `feat(make_numbat_heatmaps): add midline filtering`
- Push to origin/main automatically after meaningful commits
- Local git config uses email: `hhcg9qmx4z@privaterelay.appleid.com`

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

## Debugging Targets Pipeline Issues

When a target fails, **always parse the full dependency chain first**:
1. Identify which target(s) are errored in the log
2. Find the target definition(s) in `R/pipeline_targets_*.R` using grep: `grep -n "tar_target(TARGET_NAME" R/pipeline_targets_*.R`
3. Trace all upstream dependencies (inputs to the target function)
4. Check what each dependency **actually returns** (its return value, not just its side effects)
5. Verify downstream targets can handle what upstream targets return
6. Remember: targets tries to **serialize** the return value, not just execute the function

**Common root causes:**
- **Serialization failures**: Target returns data that can't be serialized (e.g., complex S4 objects, unevaluated expressions). Solution: return simple types (string paths, NULL, data.frame). Use `tar_file()` for file-producing targets.
- **Dependency chain breaks**: A target returns NULL or wrong type, downstream target can't handle it. Solution: trace the chain and fix the return type or have downstream read from disk directly.
- **Large data objects**: Targets tries to serialize huge data frames/lists between targets. Solution: write to disk in upstream target, have downstream read from disk instead.

**Pattern for file-producing targets:**
```r
tar_target(my_csv, {
  # do work, write file
  csv_path <- "results/output.csv"
  readr::write_csv(data, csv_path)
  # return the file path (simple string, serializes trivially)
  csv_path
})
```

For downstream targets that depend on the file:
```r
tar_target(my_plot, {
  # Depend on the file target to ensure it exists
  my_csv
  # Read from disk directly to avoid serialization issues
  data <- readr::read_csv("results/output.csv")
  plot_data(data)
})
```
