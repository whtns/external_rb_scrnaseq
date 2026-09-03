#!/bin/bash
#SBATCH --job-name=build_41_44
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_41_44_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_41_44_%j.err

# Production build for github #41, #42, #43, #44.
#
# ONE tar_make() against _targets_r431 (AGENTS.md). Checked free before
# submitting; check again before submitting anything else:
#   squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"
#
# Targets, and why these three:
#   numbat_bulk_clone_pdfs   #42/#41 -- the round-correct bulk + haplotype panels
#   sample_summaries         #42/#43 -- picks up those panels, the honest row
#                                       labels, and the new depth row
#   diploid_seu_composition  #44     -- pulls diploid_seu (the actual rebuild)
#                                       and writes the per-sample cell counts
#   numbat_rds_qc_table              -- stale (Aug 24, predates the Aug 26
#                                       majority-round change); rebuilt so the
#                                       verification table matches reality
#
# WARNING: diploid_seu OVERWRITES output/seurat/diploid_subsets/diploid_seu.rds,
# the hand-made 1.38 GB object from 2026-04-26. It was backed up first to
# diploid_seu_handmade_20260426.rds.bak (md5 8c09dbb1685a1af08b6546fc8e415de8)
# before this job was written. Every downstream diploid figure will move.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Refuse to start if the backup is not in place -- diploid_seu is destructive
# and the object it overwrites cannot be regenerated.
BAK=output/seurat/diploid_subsets/diploid_seu_handmade_20260426.rds.bak
[ -s "$BAK" ] || { echo "FATAL: backup $BAK missing; refusing to overwrite the hand-made diploid object"; exit 1; }
echo "backup present: $(ls -la $BAK)"

# install_local, not load_all: numbat_bulk_clone_pdfs, diploid_seu and
# sample_summaries all run on crew workers, which are fresh R sessions that see
# only the installed package.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers",
                                         force = TRUE, upgrade = "never", quiet = TRUE))
cat("numbatHelpers installed\n")' || { echo "FATAL: install failed"; exit 1; }

# targets does not track package function bodies. sample_summaries' COMMAND is
# unchanged, but collate_sample_summary() changed (new depth row, new labels),
# so it must be invalidated by hand or it will stay cached. Same for the qc
# table, whose inputs are unchanged but whose manifest is now stale.
# numbat_bulk_clone_pdfs and diploid_seu changed command and invalidate
# themselves, but naming them is harmless and makes the intent explicit.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(any_of(c("sample_summaries", "numbat_bulk_clone_pdfs",
                        "diploid_seu", "diploid_seu_composition",
                        "numbat_rds_qc", "numbat_rds_qc_table")))
cat("invalidated the changed targets\n")' || { echo "FATAL: invalidation failed"; exit 1; }

Rscript run_targets_bg.R \
  numbat_bulk_clone_pdfs \
  numbat_rds_qc_table \
  sample_summaries \
  diploid_seu_composition

echo "=== exit $? ==="

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"build_41_44 ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
