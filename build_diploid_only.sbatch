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

module load r/4.4.1
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
echo "### exit $?"

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"build_diploid ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
