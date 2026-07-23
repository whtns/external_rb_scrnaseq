# Pinning the numbat diploid baseline: does it stabilise 6p gain calls?

**Date:** 2026-07-22
**Jobs:** `10496424` (multi-allelic arm), `10497400` (diploid-pinned arm)
**Samples:** `SRX10264524`, `SRX14116944`
**Outputs:** `output/numbat_multiallelic/`, `output/numbat_diploidpin/`
(neither touches `output/numbat_sridhar/`, the targets pipeline's input)

## Summary

Two changes were tested against the original runs. Only one of them mattered.

- **`multi_allelic = TRUE` recovers the lost 6p gains.** `SRX10264524` lost every
  chr6p and chr1q call at round 4 of the original run; with `multi_allelic = TRUE`
  it holds four chr6p segments at every round from 2 onward, including a 21.2–41.1 Mb
  amp at LLR ~300.
- **Pinning the diploid baseline adds essentially nothing.** It reproduces the
  multi-allelic segments at the same coordinates and near-identical LLRs
  (300.8 vs 300.6), and costs a weaker chr1q call (LLR 90.0 vs 117.6).
- **The drifting-baseline hypothesis is disconfirmed.** Freezing the baseline did
  not produce convergence and did not reduce round-to-round wander. The
  segmentation itself is what is unstable, not the reference it is measured against.

The practical recommendation is therefore: **turn on `multi_allelic`, leave
`diploid_chroms = NULL`.** Do not adopt the pin.

## Background: why the pin was worth testing

`numbat:::analyze_bulk()` selects its diploid baseline through an if/else chain
(body lines 11–24):

```
exp_only | allele_only  ->  ...
diploid_chroms given    ->  use them verbatim
otherwise               ->  find_common_diploid()
```

With `diploid_chroms = NULL` the third branch runs on **every iteration**, so the
baseline against which expression fold-change is measured is re-derived each round
from the previous round's segmentation. That is a feedback loop, and the original
runs showed exactly the symptom you would predict from one: `SRX10264524` produced
12 distinct diploid sets in a single run and lost a chr1q amp at LLR 92 between
rounds 3 and 4.

Supplying `diploid_chroms` takes the second branch and `find_common_diploid()` is
never reached, making the baseline fixed and deterministic. (Note that
`common_diploid` also appears in `analyze_bulk`'s formals but appears **nowhere in
its body** — it is inert. `diploid_chroms` is the only working lever.)

Per-sample diploid sets were derived by `src/derive_diploid_chroms.R` from each
sample's own call burden across **all** consensus rounds (not just the final one),
under a hard 5% burden cap with no relaxation path, a cohort-frequency ceiling, a
config-driven veto from `large_clone_comparisons.yaml`, and chr2 banned outright.
Result: `results/diploid_chroms_per_sample.tsv`, 71 samples, 1 opted out.

- `SRX10264524` -> `3,4,5,7,8,9,11,12,14,15,18,19,22` (13 chroms, 7451 genes)
- `SRX14116944` -> `3,4,7,8,9,10,12,14,15,18,19,20,21,22` (14 chroms, 7110 genes)

## The pin verifiably took effect

| check | pinned arm | unpinned arm |
|---|---|---|
| `Using diploid chromosomes given` | 65× / 72× | 0× |
| `diploid regions:` (re-detection ran) | **0×** | 19 / 25 distinct sets |
| `No genes left in diploid regions` (silent all-gene fallback) | 0× | 0× |

So the negative results below are not a failed intervention — the intervention
did precisely what it was supposed to do.

## Result 1: `multi_allelic` recovers the 6p gains

chr6p segment calls at round 4 (the original runs' `max_iter`, so all three are
compared at equal iteration count):

**SRX10264524**

| run | chr1q | chr6p |
|---|---|---|
| original | *(none — both lost at round 4)* | *(none)* |
| multi-allelic | 1d amp LLR 117.6 (430 genes) | 6a del 6.1 · 6b amp 41.1 · 6c loh 23.3 · **6d amp 300.6** |
| diploid-pinned | 1d amp LLR 90.0 (430 genes) | 6a del 5.5 · 6b amp 41.2 · 6c loh 23.1 · **6d amp 300.8** |

The two new arms agree to within 0.1% on the chr6 21.2–41.1 Mb amp. The pin is not
doing the work; `multi_allelic` is.

**SRX14116944** never lost its chr6p call, and all three runs agree it is a large
whole-arm amp — though they disagree on its extent (0.3–49.5 Mb original,
0.3–69.8 Mb multi-allelic, 0.3–69.3 Mb pinned) and on LLR (819.3 / 801.7 / 433.9).
The pinned run's markedly lower LLR on a segment all three call is another point
against the pin.

## Result 2: the pin does not stabilise anything

Non-neutral segment count and chr6p segment count, per consensus round:

**SRX10264524**

| run | r1 | r2 | r3 | r4 | r5 | r6 | r7 | r8 |
|---|---|---|---|---|---|---|---|---|
| original — nonneu | 19 | 12 | 10 | 11 | | | | |
| original — chr6p | 3 | 2 | 2 | **0** | | | | |
| multi-allelic — nonneu | 26 | 23 | 23 | 23 | 22 | 25 | 24 | 26 |
| multi-allelic — chr6p | 3 | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
| pinned — nonneu | 26 | 23 | 23 | 25 | | | | |
| pinned — chr6p | 3 | 4 | 4 | 4 | | | | |

**SRX14116944**

| run | r1 | r2 | r3 | r4 | r5 | r6 | r7 | r8 |
|---|---|---|---|---|---|---|---|---|
| original — nonneu | 21 | 20 | 19 | 14 | | | | |
| original — chr6p | 2 | 1 | 1 | 1 | | | | |
| multi-allelic — nonneu | 31 | 31 | 26 | 23 | 28 | 28 | 28 | 30 |
| multi-allelic — chr6p | 2 | 2 | 2 | 1 | 2 | 2 | 1 | 2 |
| pinned — nonneu | 31 | 21 | 27 | 27 | | | | |
| pinned — chr6p | 2 | 1 | 2 | 1 | | | | |

Read the pinned rows against the multi-allelic rows at matched rounds. The wander
is the same size. `SRX14116944` still oscillates 2→1→2→1 on chr6p with a completely
frozen baseline, and its total call count swings 31→21→27. **Convergence was never
reached in either arm** despite `check_convergence = TRUE`; the multi-allelic run
consumed all 8 iterations without the check firing.

This is the substantive finding. The instability documented earlier — only 4 of 71
original runs settled — was attributed to the baseline feedback loop. That
attribution is wrong. Removing the loop entirely leaves the instability intact, so
its source is downstream, in the HMM segmentation and the clone-assignment step
that feeds the next round's pseudobulks.

The pin also does not pay for itself in runtime: ~60 min/iteration pinned vs
~53 min/iteration unpinned. `find_common_diploid()` is cheap next to the HMM.

## Caveat on the segment counts

`segs_consensus` can emit two rows for one interval when two consensus segments map
to it. `SRX14116944` chr1 146.9–248.9 Mb appears twice at identical LLR 752.8, once
as `amp` (`seg_cons = 1g`, `cnv_states = amp,bamp,del`, `n_states = 3`) and once as
`del` (`seg_cons = 1h`, `n_states = 1`). This is genuine multi-allelic output, not
duplication — but the tallies above count rows, so chr1q figures for that sample are
inflated by one. `n_states > 1` rows are rare: 2 in the multi-allelic round-4 file,
1 in the pinned one, 0 anywhere in the original runs.

## What to do

1. **Adopt `multi_allelic = TRUE`** for the 6p work. It is the change that recovers
   real gains, and it sharpens the calls that matter here. At round 1:
   `SRX14116944` chr1q 146.9–248.9 Mb amp **LLR 65.0 → 433.0** on identical
   coordinates; `SRX10264524` chr6p amp **66.7 → 291.3** while the segment widens
   from 21.2–33.7 Mb to 21.2–41.1 Mb (i.e. it merges 6e/6g and the intervening
   region into one stronger call rather than simply re-scoring the same interval).

   It is not a uniform gain, and should not be sold as one: `SRX14116944` chr6p
   0.3–49.5 Mb amp goes the other way, **LLR 524.9 → 465.9** on identical
   coordinates. The case for `multi_allelic` rests on recovering calls that were
   being lost outright, not on every LLR rising.
2. **Do not adopt `diploid_chroms`.** It buys no stability, weakens two calls, and
   introduces a per-sample derived parameter that has to be maintained and justified.
   `results/diploid_chroms_per_sample.tsv` and `src/derive_diploid_chroms.R` are kept
   as the record of this test, not as pipeline inputs.
3. **Iteration count remains unresolved.** `max_iter` is still a blind budget: with
   the baseline frozen the segmentation does not settle, so there is no round at
   which the answer stops changing. Any downstream use of a single round is a
   choice, not a convergence. Chasing this further means looking at the HMM
   segmentation and clone reassignment, not the reference.

## Reproducing

```bash
sbatch derive_diploid_chroms.sbatch      # writes results/diploid_chroms_per_sample.tsv
sbatch rerun_numbat_multiallelic.sbatch  # array 0-1, max_iter=8, diploid_chroms=NULL
sbatch rerun_numbat_diploidpin.sbatch    # array 0-1, max_iter=4, pinned baseline
```

Both rerun scripts differ from each other **only** in `diploid_chroms` and
`max_iter`, so the arms stay directly comparable; iteration *k* does not depend on
`max_iter`, so `segs_consensus_4.tsv` from the 4-iteration run is what the
8-iteration run would have written at round 4.

## Related

- `doc/hla_masking_6p.md` — why 6p calls are not HLA artifacts (numbat hardcodes an
  MHC blacklist in `filter_genes()`; the reference is fine)
- `src/sweep_6p_gain_candidates.R` — the missed-gain sweep, 39/39 SRX samples, no
  Tier A misses at the final consensus. Its blind spot is exactly the failure mode
  above: a gain called in an earlier round and dropped at the last is invisible to
  it, which is how `SRX10264524` passed the sweep while having lost its chr6p calls.
