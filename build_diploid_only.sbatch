#!/bin/bash
#SBATCH --job-name=build_diploid
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_diploid_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_diploid_%j.err

# github #44 only. Stage 1 (#41/#42/#43) already built and is not touched.
# ONE tar_make() against _targets_r431; checked free before submitting.
#
# WARNING: overwrites output/seurat/diploid_subsets/diploid_seu.rds (the
# hand-made 1.38 GB object). Backup guard below.

module purge
module load gcc/13.3.0   # r/4.4.1 is nested under gcc in the Lmod hierarchy;
module load r/4.4.1      # a bare purge hides it (job 11784414 failed on exactly this)
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

BAK=output/seurat/diploid_subsets/diploid_seu_handmade_20260426.rds.bak
[ -s "$BAK" ] || { echo "FATAL: backup $BAK missing"; exit 1; }
echo "backup present: $(ls -la $BAK)"

Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers",
                                         force = TRUE, upgrade = "never", quiet = TRUE))
cat("numbatHelpers installed\n")' || { echo "FATAL: install failed"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(any_of(c("diploid_seu", "diploid_seu_composition")))
cat("invalidated diploid targets\n")' || { echo "FATAL: invalidation failed"; exit 1; }

Rscript run_targets_bg.R diploid_seu_composition
echo "### run_targets_bg.R exit $?"

# run_targets_bg.R returns 0 even when tar_make() errors, so its exit status
# says nothing about whether the build worked (job 11785095 reported
# "### exit 0" on a failed diploid_seu). Ask the store directly.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
m <- tar_meta(any_of(c("diploid_seu", "diploid_seu_composition")),
              fields = c("time", "error"))
print(as.data.frame(m))
if (any(!is.na(m$error))) {
  cat("BUILD FAILED:", m$error[!is.na(m$error)], "\n"); quit(status = 1)
}
cat("diploid targets built clean\n")' || echo "### DIPLOID BUILD FAILED"

# What actually went in. A sample missing here was dropped by the
# per-sample tryCatch in assemble_diploid_seu(); grep the log for
# "diploid prep failed" to see why.
if [ -s results/diploid_seu_composition.csv ]; then
  echo '=== diploid_seu composition ==='; cat results/diploid_seu_composition.csv
fi

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"build_diploid ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
