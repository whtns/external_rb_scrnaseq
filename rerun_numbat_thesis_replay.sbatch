#!/bin/bash
#SBATCH --job-name=nb_thesis
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=0-8%3
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_thesis_%A_%a.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_thesis_%A_%a.err

# Why the three thesis-vs-current SCNA discrepancies differ, tested on the three
# tumors where they were reported. See docs/srx_vs_srr_numbat_settings.md appendix B.
#
#   SRX10264524 (RB05, was SRR13884247)  thesis: 6p gain + 11q loss
#   SRX10264525 (RB06, was SRR13884248)  thesis: 2p gain (66 Mb, LLR 222-362)
#   SRX14116947 (was SRR17960481)        thesis: 2p 85 Mb LLR 58-76 + chr9 amp 100 Mb
#
# WHAT WE ALREADY KNOW, AND WHY THIS IS NOT A THESIS-PARAMETER JOB ONLY
#
#   SRX10264524 has already been rerun twice, and BOTH reruns reproduce the
#   thesis call at PRODUCTION parameters (t=1e-2, alpha=1e-3, min_LLR=5, numbat
#   1.5.2). Comparing output/numbat_sridhar vs numbat_multiallelic vs
#   numbat_diploidpin for that sample:
#
#     thesis (SRR13884247)  6p: one amp 8.5-41.1 Mb (32 Mb) LLR 310/321
#                           11q: del 108.2-129.9 Mb (22 Mb) LLR 22/29
#
#     arm                       6p                            11q
#     production                13/8/4 Mb fragments at        absent in all
#     (multi_allelic=FALSE)     shifting coords, LLR 8-96,    4 rounds
#                               none at all by round 4
#     multi_allelic=TRUE        the 8.5-41.1 Mb span is       del 108.2-129.9 Mb,
#     (no diploid pin)          back and stable rounds 1-8,   LLR 26-32, stable
#                               but split amp 8.5-12.9 (4 Mb, rounds 1-8 --
#                               LLR 41) / loh 12.9-21.1       coordinates identical
#                               (8 Mb, LLR 24) / amp          to the thesis call
#                               21.2-41.1 (20 Mb, LLR 305)
#     multi_allelic=TRUE        same, LLR 292-305 on the      del 108.2-129.9 Mb,
#     + diploid_chroms pinned   20 Mb block, rounds 1-4       LLR 25-31, rounds 1-4
#
#   So 11q comes back exactly. 6p comes back as a stable three-part structure
#   over the thesis interval rather than the single 32 Mb amp the thesis called
#   -- a partial match, and worth saying so rather than claiming the thesis
#   segmentation was reproduced.
#
#   Both recovering arms share multi_allelic=TRUE and nothing else: the
#   multi_allelic arm has no diploid pin, so the pin is not what recovers it.
#   check_convergence differs too but only controls early stopping, not the
#   content of any round. So on this sample multi_allelic=TRUE, at unchanged
#   production parameters, restores 11q exactly and the 6p interval.
#
#   That is in tension with the rejection note in rerun_numbat_arm.sbatch
#   ("TESTED 2026-08-20 AND REJECTED ... zero canonical events gained"). That
#   verdict was scored on presence/absence of canonical events against the
#   cross-round union, which cannot see either effect at issue here:
#     - 11q_loss is not in the canonical set at all
#       (results/rb_scna_canonical_by_sample.csv carries only 13q_loss, 16q_loss,
#        1q_gain, 2p_gain, 6p_gain), so recovering it scores as nothing; and
#     - 6p_gain was already "present" in production's union as
#       recovered_marginal, so going from shifting LLR-8 fragments to a
#       stable LLR-300 block also scores as nothing.
#   The arm was rejected on a metric blind to fragmentation and to segment
#   quality. That does NOT overturn the rejection -- the over-merging seen on
#   SRX10264518 is real -- but it does mean the arm was never evaluated against
#   the question this job asks.
#
# THE THREE ARMS
#
#   ma      production parameters, multi_allelic = TRUE, nothing else changed.
#           Asks: does the SRX10264524 recovery generalise to 525 and 947?
#           SRX10264524 is included as the positive control -- it should
#           reproduce the numbat_multiallelic result above. If it does not,
#           something other than multi_allelic is in play and the other two
#           tasks are uninterpretable.
#
#   thesis  each sample's OWN thesis-era parameters, read off the workstation
#           log.txt (t=1e-5, alpha=1e-4, min_LLR=2, max_entropy=0.5, max_iter=2,
#           per-sample init_k), multi_allelic=FALSE as the thesis runs had it,
#           and the post-2026-03 cell prefilters disabled.
#           Asks: how much of the discrepancy is parameters rather than the
#           multi_allelic switch, holding the numbat version at 1.5.2?
#
#   thesis_122
#           identical to `thesis`, but run against numbat 1.2.2 -- the version
#           the workstation actually used -- instead of 1.5.2. Installed
#           2026-09-13 from the CRAN archive tarball
#           (cran.r-project.org/src/contrib/Archive/numbat/numbat_1.2.2.tar.gz,
#           4.0 MB, published 2023-02-14) into a SEPARATE library at
#           ~/R/numbat-1.2.2, reached only by setting R_LIBS in this arm.
#           The default library still holds 1.5.2 and is untouched, so the
#           production pipeline is unaffected -- verified after install.
#
#           run_numbat.R needs no modification: every argument it passes is in
#           1.2.2's run_numbat signature, and both internals it reaches for
#           (numbat:::get_bulk, numbat:::detect_clonal_loh) exist there.
#           1.2.2's defaults differ (min_cells 50, init_k 3, max_entropy 0.5,
#           multi_allelic TRUE) but the script passes all of them explicitly.
#
#           `thesis` vs `thesis_122` is the version effect with parameters held
#           fixed; `thesis_122` vs the archived workstation output is as close
#           to a like-for-like replay as this machine can get.
#
# WHAT NO ARM CONTROLS
#
#   The cell set. SRR13884247 and SRR13884248 were run on the workstation
#   against a RANDOM 10,000-cell subsample under the old cell_ceiling, which is
#   not reproducible (no seed was recorded, and cell_ceiling is commented out of
#   run_numbat.R). Setting min_allele_depth=0 / min_snps_per_cell=0 in both
#   thesis arms removes the prefilters added in commit 162eefd, which is the
#   part of the cell-set difference that IS reproducible; the subsample is not.
#   So even thesis_122 can differ from the archived workstation run on cell
#   membership alone, and a residual difference there is not evidence of a code
#   or parameter effect.
#
# WRITES
#   output/numbat_thesis_replay/<arm>/<sample>/ only.
#   output/numbat_sridhar/ (the targets pipeline's input) is NOT touched, nor are
#   the existing numbat_multiallelic/ and numbat_diploidpin/ arms.
#   None of the three samples is one of the held-out SRX11133592/93/94.

set -euo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj

# task -> (sample, arm). Grouped by arm so a partial array runs whole arms:
#   --array=0-2  ma only      --array=3-5  thesis only     --array=6-8  thesis_122 only
SAMPLES=(SRX10264524 SRX10264525 SRX14116947 \
         SRX10264524 SRX10264525 SRX14116947 \
         SRX10264524 SRX10264525 SRX14116947)
ARMS=(ma ma ma  thesis thesis thesis  thesis_122 thesis_122 thesis_122)

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
ARM=${ARMS[$SLURM_ARRAY_TASK_ID]}
OUTDIR=$PROJ/output/numbat_thesis_replay/$ARM/$S

# Guard: never write into the production store, whatever anyone edits above.
case "$OUTDIR" in
  *numbat_sridhar*) echo "FATAL: refusing to write into numbat_sridhar" >&2; exit 1 ;;
esac

# Per-sample parameters.
#   ma         = production values, verbatim from output/numbat_sridhar/<S>/log.txt
#   thesis     = workstation values, verbatim from the matched SRR log.txt
#   thesis_122 = same values as thesis; only the numbat library differs
case "$S:$ARM" in
  # production: t 1e-2, alpha 1e-3, gamma 20, init_k 5, max_iter 4, min_LLR 5, max_entropy 0.7
  SRX10264524:ma)     T=1e-2; ALPHA=1e-3; GAMMA=20; INIT_K=5; MAX_ITER=4; MIN_LLR=5; MAX_ENT=0.7; MA=TRUE;  AD=5; SNP=50 ;;
  SRX10264525:ma)     T=1e-2; ALPHA=1e-3; GAMMA=20; INIT_K=5; MAX_ITER=4; MIN_LLR=5; MAX_ENT=0.7; MA=TRUE;  AD=5; SNP=50 ;;
  SRX14116947:ma)     T=1e-2; ALPHA=1e-3; GAMMA=20; INIT_K=5; MAX_ITER=4; MIN_LLR=5; MAX_ENT=0.7; MA=TRUE;  AD=5; SNP=50 ;;
  # thesis: from SRR13884247 / SRR13884248 / SRR17960481 log.txt (numbat 1.2.2).
  # alpha was never passed at the old call site, so 1e-4 is numbat's default, not a choice.
  # init_k differs per sample and no longer matches metadata.tsv (num_clone = 5 for all three).
  SRX10264524:thesis) T=1e-5; ALPHA=1e-4; GAMMA=20; INIT_K=6; MAX_ITER=2; MIN_LLR=2; MAX_ENT=0.5; MA=FALSE; AD=0; SNP=0 ;;
  SRX10264525:thesis) T=1e-5; ALPHA=1e-4; GAMMA=20; INIT_K=3; MAX_ITER=2; MIN_LLR=2; MAX_ENT=0.5; MA=FALSE; AD=0; SNP=0 ;;
  SRX14116947:thesis) T=1e-5; ALPHA=1e-4; GAMMA=20; INIT_K=5; MAX_ITER=2; MIN_LLR=2; MAX_ENT=0.5; MA=FALSE; AD=0; SNP=0 ;;
  # thesis_122: identical parameters, numbat 1.2.2 instead of 1.5.2
  SRX10264524:thesis_122) T=1e-5; ALPHA=1e-4; GAMMA=20; INIT_K=6; MAX_ITER=2; MIN_LLR=2; MAX_ENT=0.5; MA=FALSE; AD=0; SNP=0 ;;
  SRX10264525:thesis_122) T=1e-5; ALPHA=1e-4; GAMMA=20; INIT_K=3; MAX_ITER=2; MIN_LLR=2; MAX_ENT=0.5; MA=FALSE; AD=0; SNP=0 ;;
  SRX14116947:thesis_122) T=1e-5; ALPHA=1e-4; GAMMA=20; INIT_K=5; MAX_ITER=2; MIN_LLR=2; MAX_ENT=0.5; MA=FALSE; AD=0; SNP=0 ;;
  *) echo "FATAL: no parameter row for $S:$ARM" >&2; exit 1 ;;
esac

mkdir -p "$OUTDIR"

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

# numbat version selection. The default library holds 1.5.2; the 1.2.2 build
# lives in its own tree and is reached ONLY by prepending it to R_LIBS here.
# The default library stays on the path so Seurat, data.table, fs, dplyr etc.
# still resolve -- only numbat is shadowed.
DEFAULT_RLIB=$HOME/R/x86_64-pc-linux-gnu-library/4.4
NUMBAT_122_LIB=$HOME/R/numbat-1.2.2
if [ "$ARM" = "thesis_122" ]; then
  if [ ! -d "$NUMBAT_122_LIB/numbat" ]; then
    echo "FATAL: numbat 1.2.2 not installed at $NUMBAT_122_LIB" >&2
    echo "       install: Rscript -e \"install.packages('numbat_1.2.2.tar.gz', lib='$NUMBAT_122_LIB', repos=NULL, type='source')\"" >&2
    exit 1
  fi
  export R_LIBS="$NUMBAT_122_LIB:$DEFAULT_RLIB"
else
  export R_LIBS="$DEFAULT_RLIB"
fi

cd $PROJ

echo "=== $S arm=$ARM -> $OUTDIR ==="
echo "t=$T alpha=$ALPHA gamma=$GAMMA init_k=$INIT_K max_iter=$MAX_ITER min_LLR=$MIN_LLR"
echo "max_entropy=$MAX_ENT multi_allelic=$MA min_allele_depth=$AD min_snps_per_cell=$SNP"
# Assert the intended numbat is the one that will load -- a silent fallback to
# 1.5.2 would make the thesis_122 arm a duplicate of thesis without saying so.
NB_VER=$(Rscript -e 'cat(as.character(packageVersion("numbat")))' 2>/dev/null)
NB_LOC=$(Rscript -e 'cat(dirname(find.package("numbat")))' 2>/dev/null)
echo "numbat: $NB_VER  from $NB_LOC"
if [ "$ARM" = "thesis_122" ] && [ "$NB_VER" != "1.2.2" ]; then
  echo "FATAL: arm thesis_122 resolved numbat $NB_VER, expected 1.2.2" >&2; exit 1
fi
if [ "$ARM" != "thesis_122" ] && [ "$NB_VER" != "1.5.2" ]; then
  echo "FATAL: arm $ARM resolved numbat $NB_VER, expected 1.5.2" >&2; exit 1
fi
echo "start: $(date)"

# Not fatal under set -e: a numbat failure here would skip the readout below and
# leave a multi-hour run with no diagnostics. Capture, report, exit with it.
set +e
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
  "max_entropy=$MAX_ENT" \
  "min_LLR=$MIN_LLR" \
  "max_iter=$MAX_ITER" \
  "tau='0.3'" \
  "read_prop=1" \
  "multi_allelic=$MA" \
  "check_convergence=FALSE" \
  "min_allele_depth=$AD" \
  "min_snps_per_cell=$SNP" \
  "ncores=8"
RC=$?

echo "numbat exit: $RC ; end: $(date)"

# Readout: did the thesis events come back? Print every non-neutral segment on
# the chromosomes at issue, per round, so the answer is in the job log.
#   6p  chr6  < 59.8 Mb      11q chr11 > 53.4 Mb
#   2p  chr2  < 93.9 Mb      chr9 whole
echo "=== thesis events, all rounds ==="
python3 - "$OUTDIR" "$S" <<'PY'
import csv, glob, os, sys
outdir, samp = sys.argv[1], sys.argv[2]
WANT = {"SRX10264524": [("6p",6,0,59.8), ("11q",11,53.4,135.1)],
        "SRX10264525": [("2p",2,0,93.9)],
        "SRX14116947": [("2p",2,0,93.9), ("chr9",9,0,138.4)]}
for f in sorted(glob.glob(os.path.join(outdir, "segs_consensus_*.tsv"))):
    hits = []
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["cnv_state_post"] in ("neu", ""):
            continue
        ch = int(r["CHROM"]); s = float(r["seg_start"])/1e6; e = float(r["seg_end"])/1e6
        for lbl, c, lo, hi in WANT.get(samp, []):
            # any overlap with the arm window, segment at least 1 Mb
            if ch == c and e > lo and s < hi and (e - s) >= 1.0:
                hits.append(f"{lbl} {r['cnv_state_post']} {s:.1f}-{e:.1f}Mb "
                            f"({e-s:.0f}Mb) LLR={float(r['LLR']):.0f}")
    print(f"-- {os.path.basename(f)}: " + ("; ".join(hits) if hits else "NONE"))
PY

echo "=== parameters numbat actually used ==="
sed -n '/Running under parameters/,/Input metrics/p' "$OUTDIR/log.txt" 2>/dev/null || \
  echo "no log.txt at $OUTDIR"

exit $RC
