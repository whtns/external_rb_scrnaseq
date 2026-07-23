#!/bin/bash
#SBATCH --job-name=nb_diploidpin
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_diploidpin_%A_%a.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_diploidpin_%A_%a.err

# Pin the diploid baseline per sample and rerun, to test whether the iteration
# instability is caused by the baseline drifting.
#
# Motivation: analyze_bulk's diploid logic is an if/else chain (body lines 11-24).
# With diploid_chroms = NULL it falls through to find_common_diploid(), which
# re-derives the baseline EVERY iteration. SRX10264524 produced 12 distinct
# diploid sets in one run and lost a chr1q amp at LLR 92 between rounds 3 and 4.
# Supplying diploid_chroms takes the earlier branch and never reaches
# find_common_diploid(), making the baseline fixed and deterministic.
#
# Differs from rerun_numbat_multiallelic.sbatch (job 10496424) ONLY by
# diploid_chroms and max_iter, so the two remain directly comparable: both use
# multi_allelic=TRUE, check_convergence=TRUE, and otherwise the original
# per-sample parameters.
#
# max_iter = 4 to allow evaluation sooner, and it costs nothing in comparability:
# iteration k does not depend on max_iter, so segs_consensus_4.tsv here is
# identical to what an 8-iteration run would have written at round 4. It also
# matches the ORIGINAL runs' max_iter, so the three-way comparison below is at
# equal iteration count rather than converged-vs-unconverged. If
# check_convergence fires before round 4, the run stops earlier still.
#
# Writes to output/numbat_diploidpin/<sample>/ -- neither output/numbat_sridhar/
# (the targets pipeline's input) nor output/numbat_multiallelic/ is touched.

set -euo pipefail

SAMPLES=(SRX10264524 SRX14116944)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
TSV=$PROJ/results/diploid_chroms_per_sample.tsv
OUTDIR=$PROJ/output/numbat_diploidpin/$S

if [ ! -f "$TSV" ]; then
  echo "FATAL: $TSV missing -- run 'sbatch derive_diploid_chroms.sbatch' first" >&2
  exit 1
fi

# Column 2 of the row for this sample. "NA" means the sample opted out: it could
# not assemble an adequate CLEAN baseline under the hard burden cap, so pinning it
# would mean designating CNV-laden chromosomes as diploid. Refuse rather than run
# a knowably wrong baseline.
DIP=$(awk -F'\t' -v s="$S" 'NR>1 && $1==s {print $2}' "$TSV")
if [ -z "$DIP" ]; then
  echo "FATAL: no diploid_chroms row for $S in $TSV" >&2
  exit 1
fi
if [ "$DIP" = "NA" ]; then
  echo "FATAL: $S opted out -- no adequate clean baseline under the burden cap." >&2
  echo "       Pinning it would designate CNV-laden chromosomes as diploid." >&2
  echo "       Leave this sample on numbat's default find_common_diploid()." >&2
  exit 1
fi

mkdir -p "$OUTDIR"

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

cd $PROJ

echo "=== $S : diploid_chroms = $DIP -> $OUTDIR ==="
echo "start: $(date)"

# Deliberately not fatal: under `set -e` a numbat failure would abort the script
# before the verification block below, so a multi-hour run could die and leave no
# diagnostics. Capture the status, always run the checks, exit with it at the end.
set +e
Rscript pipeline/scripts/run_numbat.R \
  "matrix_file='$PROJ/output/cellranger/$S/outs/filtered_feature_bc_matrix/matrix.mtx.gz'" \
  "seu_path='$PROJ/output/seurat/${S}_seu.rds'" \
  "allele_df='$PROJ/output/numbat/${S}_allele_counts.tsv.gz'" \
  "ref_path='/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds'" \
  "out_dir='$OUTDIR'" \
  "init_k=5" \
  "t='1e-2'" \
  "alpha='1e-3'" \
  "gamma=20" \
  "max_entropy=0.7" \
  "min_LLR=5" \
  "max_iter=4" \
  "tau='0.3'" \
  "read_prop=1" \
  "multi_allelic=TRUE" \
  "check_convergence=TRUE" \
  "diploid_chroms='$DIP'" \
  "ncores=8"
RC=$?

echo "numbat exit: $RC ; end: $(date)"

# Did the pin actually take effect, and did the baseline survive?
echo "=== baseline sanity ==="
LOG="$OUTDIR/log.txt"
if [ -f "$LOG" ]; then
  echo "-- 'Using diploid chromosomes given' (expect one per analyze_bulk call) --"
  grep -c "Using diploid chromosomes given" "$LOG" || echo "  0 -- PIN DID NOT TAKE EFFECT"
  grep -m1 "Using diploid chromosomes given" "$LOG" || true
  echo "-- find_common_diploid should NEVER appear --"
  grep -c "diploid regions:" "$LOG" || echo "  0 (correct: re-detection did not run)"
  echo "-- empty-baseline fallback (silently uses ALL genes; must be absent) --"
  grep -c "No genes left in diploid regions" "$LOG" || echo "  0 (correct)"
  echo "-- convergence --"
  grep -E "Iteration|Convergence reached" "$LOG" | tail -12 || true
fi

# chr6p and chr1 calls per consensus round: original vs multi_allelic vs pinned.
echo "=== $S calls by consensus round, three runs ==="
for D in "$PROJ/output/numbat_sridhar/$S" "$PROJ/output/numbat_multiallelic/$S" "$OUTDIR"; do
  echo "--- $D ---"
  for f in "$D"/segs_consensus_*.tsv; do
    [ -f "$f" ] || continue
    awk -F'\t' -v FN="$(basename "$f")" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
      $h["cnv_state"]!="neu" && $h["cnv_state"]!="" {
        n++;
        if($h["CHROM"]==6 && $h["seg_start"]+0<58500000) p6++;
        if($h["CHROM"]==1 && $h["seg_start"]+0>125000000) q1++
      } END{printf "  %-22s non-neutral=%3d  chr6p=%d  chr1q=%d\n", FN, n+0, p6+0, q1+0}' "$f"
  done
done

exit $RC
