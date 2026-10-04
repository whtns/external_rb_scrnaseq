#!/bin/bash
#SBATCH --job-name=nb_merged_reps
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=240G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_merged_reps_%A_%a.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_merged_reps_%A_%a.err

# numbat on the Wu et al. technical-replicate pairs, merged per tumor.
# Sensitivity check ALONGSIDE the cohort: writes only output/numbat_merged/<tumor>/.
# Nothing under output/numbat/, numbat_sridhar/ or numbat_t1e5/ is touched, and
# the 6 replicate samples stay in the targets cohort as they are.
#
#   RB01 = SRX10264517 + SRX10264518
#   RB02 = SRX10264519 + SRX10264520
#   RB03 = SRX10264521 + SRX10264522
#
# The captures are independent (~1% barcode overlap = whitelist collisions), so
# cells are concatenated with "<SRX>_" barcode prefixes. Stages:
#   1. joint genotype + Eagle phasing from the existing per-replicate pileups
#      (container: numbat 1.5.2 + eagle, the stack that made the originals)
#   2. merged count matrix / metadata / allele table
#   3. run_numbat.R with the production t=1e-5 parameters (identical in all six
#      replicates' log.txt: t=1e-5 alpha=1e-3 gamma=20 init_k=5 min_LLR=5
#      max_entropy=0.7 tau=0.3 multi_allelic=FALSE, max_iter=4)
# Each stage is skipped if its output already exists, so a rerun resumes.
#
# 240G: single captures ran in 128G at <=7.5k cells; merged tumors are 6-15k.

set -euo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
SIF=$PROJ/pipeline/containers/numbat-pipeline.sif
cd "$PROJ"

TUMORS=(RB01 RB02 RB03)
PAIRS=(SRX10264517,SRX10264518 SRX10264519,SRX10264520 SRX10264521,SRX10264522)
TUMOR=${TUMORS[$SLURM_ARRAY_TASK_ID]}
PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
OUT=$PROJ/output/numbat_merged/$TUMOR
mkdir -p "$OUT"

export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

echo "=== $TUMOR ($PAIR) -> $OUT  start: $(date)"

# 1. joint phasing
if [ ! -f "$OUT/${PAIR##*,}_allele_counts.tsv.gz" ]; then
  module purge
  module load apptainer
  apptainer exec "$SIF" \
    Rscript src/phase_merged_replicates.R "${TUMOR}_merged" "$PAIR" "$OUT" 8
fi

module purge
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# 2. merged inputs
if [ ! -f "$OUT/allele_counts.tsv.gz" ]; then
  Rscript src/build_merged_replicate_inputs.R "$PAIR" "$OUT"
fi

# 3. numbat
NB=$OUT/numbat
mkdir -p "$NB"
echo "numbat start: $(date)"
Rscript pipeline/scripts/run_numbat.R \
  "matrix_file='$OUT/matrix/matrix.mtx.gz'" \
  "seu_path='$OUT/seu_meta.rds'" \
  "allele_df='$OUT/allele_counts.tsv.gz'" \
  "ref_path='/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds'" \
  "out_dir='$NB'" \
  "init_k=5" \
  "t='1e-5'" \
  "alpha='0.001'" \
  "gamma=20" \
  "max_entropy=0.7" \
  "min_LLR=5" \
  "max_iter=4" \
  "tau='0.3'" \
  "read_prop=1" \
  "multi_allelic=FALSE" \
  "check_convergence=FALSE" \
  "ncores=8"
echo "end: $(date)"
