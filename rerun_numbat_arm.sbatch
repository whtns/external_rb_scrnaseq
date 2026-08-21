#!/bin/bash
#SBATCH --job-name=nb_arm
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=0-35%8
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_arm_%A_%a.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_arm_%A_%a.err

# Rerun numbat over the rerunnable cohort with multi_allelic = TRUE, optionally
# overriding t. Generalises rerun_numbat_multiallelic.sbatch from 2 samples to
# the whole SRX set, and replays each sample's OWN logged parameters so the only
# differences from its production run are the ones named below.
#
# WHY EACH SAMPLE'S PARAMS COME FROM ITS log.txt, NOT config.yaml:
#   pipeline/config.yaml has drifted from what production actually ran.
#   config says alpha 1e-4; 31 of 39 SRX samples ran at alpha 1e-3. numbat echoes
#   the real values under "Running under parameters:" in
#   output/numbat_sridhar/<sample>/log.txt, and those are authoritative.
#   See results/numbat_param_regimes.csv (src/diag_numbat_param_regimes.R).
#
# WHY ONLY 39 SAMPLES:
#   The 32 SRR samples were run on the previous machine and their inputs
#   (cellranger matrix, seurat object, allele counts) were never migrated here,
#   so they cannot be rerun without regenerating from FASTQ. Only the 39 SRX
#   samples have all three inputs on disk. The cross-round union in
#   src/diag_rb_scna_recovery.R is the only recovery lever that reaches the SRR
#   samples at all.
#
# DELIBERATE CHANGES vs each sample's production run:
#   multi_allelic  FALSE -> TRUE   numbat's own default; production hardcoded
#                                  FALSE at the run_numbat() call site.
#
#                                  !! TESTED 2026-08-20 AND REJECTED. Do not
#                                  promote this arm. Pilot job 11228996, 5
#                                  scorable samples: union 9 vs production's 10.
#                                  Zero canonical events gained, one lost
#                                  (SRX10264518 6p). TRUE over-merges -- chr6
#                                  collapses to one 65-170Mb loh/del block with
#                                  runaway LLR (288->429->1513->1585) that
#                                  swallows the focal 6p amp.
#
#                                  Two earlier claims here were wrong:
#                                  (a) FALSE does NOT disable test_multi_allelic()
#                                      -- 24/71 production dirs already have
#                                      bamp/bdel under FALSE;
#                                  (b) the "SRX10264524 1 -> 4, SRX14116944 3 -> 4"
#                                      validation compared against production's
#                                      FINAL ROUND, not its cross-round union.
#                                      Against the union it is 5 -> 4 and 4 -> 4.
#                                  docs/numbat_scna_recovery_recommendations.md sec 3a
#   max_iter       -> 4            37/39 SRX samples already ran at 4; this
#                                  pins the 2 stragglers to match. Do NOT raise
#                                  it: rounds 5-8 were tested and bought only
#                                  wander at ~3.5 h/sample.
#   t              -> $T_OVERRIDE  unset by default (replays the logged t).
#
# NOT changed, all previously tested and settled:
#   diploid_chroms     stays NULL   pinning was negative; doc/diploid_baseline_pin_6p.md
#   check_convergence  stays FALSE  never fires; only 4/71 runs reach d_k = 0
#   HLA/MHC masking    n/a          numbat hardcodes it in filter_genes()
#
# USAGE
#   # Arm A -- multi_allelic only, logged t:
#   sbatch --export=ALL,ARM_NAME=ma rerun_numbat_arm.sbatch
#
#   # Arm B -- multi_allelic + numbat's default t, on a 6-sample pilot:
#   sbatch --array=0-5 --export=ALL,ARM_NAME=ma_t1e5,T_OVERRIDE=1e-5 \
#          rerun_numbat_arm.sbatch
#
#   Output goes to output/numbat_<ARM_NAME>/<sample>/. output/numbat_sridhar/
#   is the targets pipeline's input and is NEVER touched.

set -uo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
ARM_NAME=${ARM_NAME:-ma}
T_OVERRIDE=${T_OVERRIDE:-}
SRC_DIR=$PROJ/output/numbat_sridhar

# Samples held out of the arm at the user's request. Kept as an explicit list so
# the array index -> sample mapping is reproducible from the script alone.
EXCLUDE="SRX11133592 SRX11133593 SRX11133594"

cd "$PROJ"

# Rerunnable samples: an SRX dir with a production log.txt and all three inputs.
mapfile -t SAMPLES < <(
  for d in "$SRC_DIR"/SRX*/; do
    s=$(basename "$d")
    [ -f "$d/log.txt" ] || continue
    [ -f "$PROJ/output/cellranger/$s/outs/filtered_feature_bc_matrix/matrix.mtx.gz" ] || continue
    [ -f "$PROJ/output/seurat/${s}_seu.rds" ] || continue
    [ -f "$PROJ/output/numbat/${s}_allele_counts.tsv.gz" ] || continue
    case " $EXCLUDE " in *" $s "*) continue ;; esac
    echo "$s"
  done | sort
)

echo "rerunnable samples: ${#SAMPLES[@]}"
if [ "${SLURM_ARRAY_TASK_ID:-0}" -ge "${#SAMPLES[@]}" ]; then
  echo "array index ${SLURM_ARRAY_TASK_ID} beyond sample list; nothing to do"
  exit 0
fi
S=${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}

# Pull this sample's real parameters out of numbat's own echo block.
logparam () {
  awk -v k="$1" '
    /Running under parameters:/ {inblock=1; next}
    /Input metrics:/ {inblock=0}
    inblock && $1 == k && $2 == "=" {print $3; exit}
  ' "$SRC_DIR/$S/log.txt"
}

T=$(logparam t)
ALPHA=$(logparam alpha)
GAMMA=$(logparam gamma)
INIT_K=$(logparam init_k)
MIN_LLR=$(logparam min_LLR)
MAX_ENTROPY=$(logparam max_entropy)
TAU=$(logparam tau)

# init_k is echoed as a float ("5.0"); run_numbat wants an integer.
INIT_K=${INIT_K%.*}
[ -n "$T_OVERRIDE" ] && T=$T_OVERRIDE

OUTDIR=$PROJ/output/numbat_${ARM_NAME}/$S
mkdir -p "$OUTDIR"

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

echo "=== $S -> $OUTDIR (arm: $ARM_NAME) ==="
echo "replayed from log.txt: t=$T alpha=$ALPHA gamma=$GAMMA init_k=$INIT_K"
echo "                       min_LLR=$MIN_LLR max_entropy=$MAX_ENTROPY tau=$TAU"
echo "changed:               multi_allelic=TRUE max_iter=4${T_OVERRIDE:+ t=$T_OVERRIDE}"
echo "start: $(date)"

Rscript pipeline/scripts/run_numbat.R \
  "matrix_file='$PROJ/output/cellranger/$S/outs/filtered_feature_bc_matrix/matrix.mtx.gz'" \
  "seu_path='$PROJ/output/seurat/${S}_seu.rds'" \
  "allele_df='$PROJ/output/numbat/${S}_allele_counts.tsv.gz'" \
  "ref_path='/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds'" \
  "out_dir='$OUTDIR'" \
  "init_k=$INIT_K" \
  "t='$T'" \
  "alpha='$ALPHA'" \
  "gamma=$GAMMA" \
  "max_entropy=$MAX_ENTROPY" \
  "min_LLR=$MIN_LLR" \
  "max_iter=4" \
  "tau='$TAU'" \
  "read_prop=1" \
  "multi_allelic=TRUE" \
  "check_convergence=FALSE" \
  "ncores=8"

RC=$?
echo "exit: $RC ; end: $(date)"

# Immediate readout: canonical RB events, production vs this arm, per round.
echo "=== $S canonical RB events by consensus round ==="
for D in "$SRC_DIR/$S" "$OUTDIR"; do
  echo "--- $D ---"
  for f in "$D"/segs_consensus_*.tsv; do
    [ -f "$f" ] || continue
    ev=$(awk -F'\t' '
      NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
      {
        ch=$h["CHROM"]; st=$h["seg_start"]+0; en=$h["seg_end"]+0;
        s=(("cnv_state_post" in h) && $h["cnv_state_post"]!="") ? $h["cnv_state_post"] : $h["cnv_state"];
        g=(s ~ /amp/); l=(s ~ /del|loh/); e="";
        if      (ch=="1"  && en>125000000 && g) e="1q+";
        else if (ch=="2"  && st<93000000  && g) e="2p+";
        else if (ch=="6"  && st<59000000  && g) e="6p+";
        else if (ch=="13" && l)                 e="13q-";
        else if (ch=="16" && en>36800000 && l)  e="16q-";
        if (e!="") seen[e]=1;
      }
      END {n=0; o=""; for(k in seen){n++; o=o k " "} printf "%d  %s", n, o}
    ' "$f")
    echo "  $(basename "$f"): $ev"
  done
done

exit $RC
