# Choosing the number of numbat iterations: evidence for a fixed-round policy

**Scope.** 36 SRX samples (`SRX11133592/93/94` held out per project policy).
33 have four consensus rounds on disk, 3 have two. Canonical RB events use the
arm-coverage-thresholded definition (`RB_MIN_ARM_FRAC = 0.15`) that the reported
QC table `numbat_rds_qc.csv` uses.

**Bottom line.** Two of the three premises hold. The third does not:
**numbat does not converge within four rounds** (§3), so no round can be
defended as the converged answer.

A fixed default of `i = 2` was the initial recommendation, but §4b supersedes it:
round 2 is itself the sole dissenter from all other rounds in 6 sample x event
cases, and contributes only 3 unique calls of its own. **A majority-of-rounds
criterion is preferred** -- it is count-neutral against `i = 2` (76 vs 77 events),
selects a better-supported set, and sidesteps the convergence problem entirely,
since "most iterations call it" is a stability claim rather than a fixed-point
claim. Its open cost is 17 round-1-only calls it would reject, which need the
per-cell posterior check to adjudicate (§4c).

---

## 1. Round 1 over-calls segments — confirmed

| round | samples | mean segs | median segs | mean focal (<10 Mb) | median seg size |
|-------|---------|-----------|-------------|---------------------|-----------------|
| 1 | 36 | 23.8 | 21.5 | 12.9 | 7.1 Mb |
| 2 | 36 | 18.9 | 16.5 |  9.5 | 8.2 Mb |
| 3 | 33 | 20.3 | 14.0 |  9.9 | 10.2 Mb |
| 4 | 33 | 19.3 | 13.0 |  9.3 | 9.0 Mb |

Round 1 calls **858** non-neutral segments against **681** at round 2 — 26% more
— and does so in 30 of 36 samples (Wilcoxon signed-rank, *p* = 3.5 × 10⁻⁶; 3 tied,
3 reversed). The excess is concentrated in short segments: 466 focal calls at
round 1 versus 343 at round 2, and median segment size *rises* from 7.1 to
8.2 Mb as the extra calls drop out.

→ `results/round_convergence_by_round.csv`, `results/round1_vs_round2_segments.csv`

## 2. The segments round 1 adds are the low-confidence ones — confirmed

Only **408 of 858 (48%)** round-1 segments survive to round 2 with ≥50%
reciprocal overlap in the same direction. The 450 that do not are systematically
weaker:

| round-1 segments | n | median size | median LLR | median genes | % focal |
|------------------|---|-------------|------------|--------------|---------|
| dropped at r2 | 450 | 6.04 Mb | 8.3 | 31.0 | 61% |
| kept at r2 | 408 | 11.63 Mb | 9.9 | 47.5 | 47% |

Size *p* = 1.9 × 10⁻⁴, LLR *p* = 1.3 × 10⁻¹⁰ (Wilcoxon rank-sum). This is the
quantitative basis for treating round 1's surplus as enriched for artefact:
it is smaller, gene-poorer, and lower-likelihood than the segments that persist.

→ `results/round1_segment_persistence.csv`

## 3. There is no convergence point — **not** confirmed

Agreement between consecutive rounds, measured as the Jaccard index over
non-neutral segment sets (≥50% reciprocal overlap, same chromosome, same
direction):

| pair | samples | mean Jaccard | mean kept | mean new |
|------|---------|--------------|-----------|----------|
| r1→r2 | 36 | 0.355 | 0.447 | 0.376 |
| r2→r3 | 33 | 0.490 | 0.636 | 0.356 |
| r3→r4 | 33 | 0.571 | 0.685 | 0.259 |

Agreement rises monotonically and is **still rising** at the last available step
(r2→r3 vs r3→r4, *p* = 0.039). Only **4 of 33** samples have an unchanged segment
set between rounds 3 and 4; **zero** are stable at r1→r2 or r2→r3. Even at the
final step, roughly a quarter of segments are new and a third of the previous
round's segments have been dropped.

**Implication for the manuscript.** Do not write that numbat has converged, and
do not justify the final round on those grounds either — the final round is no
more settled than round 2. The honest framing is that numbat's iteration does
not reach a fixed point on this data within four rounds.

Two conclusions follow from that, and only the first survives §4b. Either the
reported round is *chosen a priori and applied uniformly* — which would favour
numbat's own documented default `Numbat$new(..., i = 2)` as the choice requiring
least special pleading — or the reporting criterion is made independent of any
single round. §4b shows the first option fails on its own terms, because round 2
is not a neutral reference: it is the sole dissenter in 6 cases. The second
option, a majority-of-rounds criterion, is what §4c recommends, and it turns the
absence of convergence from a weakness into a non-issue: reporting an event when
most iterations call it is a stability statement that never needs a fixed point.

→ `results/round_convergence_stability.csv`

## 4. What a fixed `i = 2` costs

| | events |
|---|---|
| canonical RB events at `i = 2` | 77 |
| canonical RB events, best round per sample | 104 |
| missed by `i = 2` | **27 (26%)** |
| samples affected | **20 of 36** |

Missed by event: 16q loss (7), 6p gain (6), 13q loss (5), 2p gain (5), 1q gain (4).

**All 27 are recovered at round 1**; 17 of the 27 at round 1 *only*. The other 10
are also callable at round 3 or 4. The deviation rule is therefore single-valued
in practice and needs no search: *if a canonical RB locus is not called at
`i = 2`, check `i = 1`* — round 1 recovers every one of them, and no deeper round
does.

### Why not iterate further instead?

Going past round 2 does not accumulate RB events — it churns them:

| round | samples | RB events | gained vs r2 | lost vs r2 | net |
|---|---|---|---|---|---|
| 1 | 36 | 98 | 27 (20 samples) | 6 (5 samples) | **+21** |
| 2 | 36 | 77 | — | — | — |
| 3 | 33 | 76 | 9 (9 samples) | 10 (7 samples) | −1 |
| 4 | 33 | 71 | 7 (7 samples) | 13 (10 samples) | −6 |

Rounds 3 and 4 call *fewer* canonical events than round 2 in total, and their
sets are not nested with it: each gains a handful the previous round lacked while
dropping a comparable or larger number. Rounds 3 and 4 combined recover only 10
of the 27 events missed at `i = 2` (5 × 16q loss, 4 × 6p gain, 1 × 1q gain).
Round 1 is the only round that adds canonical events on net, which is the second
reason the deviation points there rather than deeper.

Splitting the 27 by what round 2 actually shows at the locus:

| status at r2 | n | median arm coverage at r2 | median arm coverage at r1 | median LLR at r1 |
|---|---|---|---|---|
| absent (no call on the arm) | 17 | 0.000 | 0.701 | 94.7 |
| sub-threshold (called, but <15% of arm) | 10 | 0.098 | 0.530 | 75.3 |

This is the key point for defending the deviation. The round-1 calls being
recovered are **not** the focal, low-LLR class from §2 — those had median size
6 Mb and median LLR 8.3. The recovered RB events cover a median **53–70% of the
target arm** with median LLR **75–95**. They are arm-scale, high-confidence
calls. Deviating to round 1 at a named RB locus does not import the artefact
population that makes round 1 unattractive as a global default.

Ten of the 27 are not disagreements about whether the arm is altered at all —
round 2 calls the arm, but fragmented below the 15% coverage floor.

→ `results/round_policy_missed_events.csv`, `results/round_policy_i2_vs_union.csv`

## 4b. Round 2 is itself sometimes the outlier

Sections 1–4 all use round 2 as the reference for every contrast, which by
construction cannot reveal round 2 being the round in error. Reading across
rounds instead shows that it sometimes is.

**Round 2 alone dissents from all other rounds in 6 sample × event cases:**

| sample | event | rounds calling | mean arm coverage where called |
|---|---|---|---|
| SRX10264522 | 16q_loss | 1,3,4 | 0.718 |
| SRX10264526 | 16q_loss | 1,3,4 | 0.265 |
| SRX10831282 | 16q_loss | 1,3,4 | 0.458 |
| SRX14116946 | 6p_gain | 1,3,4 | 0.653 |
| SRX22868103 | 6p_gain | 1,3,4 | 0.471 |
| SRX22868105 | 6p_gain | 1,3,4 | 0.495 |

Ten cases in total have ≥2 other rounds calling an event round 2 misses.
SRX10264526 is the clearest: rounds 1, 3 and 4 all call essentially the same
chr16 interval (~69.8–81.1 Mb, ~11 Mb, ~40 genes, del, LLR 7.8–8.5), and at
round 2 there is **no del/loh segment anywhere on 16q** — coverage 0.000, not a
sub-threshold fragment.

Conversely, round 2 uniquely adds only **3** events (SRX11133585 13q, SRX14116947
16q, SRX22868103 13q), two of them at weak arm coverage (0.16, 0.41) — the
profile of a round-specific artefact rather than a recovered event.

→ `results/rb_event_by_round_matrix.csv`

## 4c. Policy comparison

| policy | canonical RB events |
|---|---|
| fixed `i = 2` | 77 |
| **majority of available rounds** | **76** |
| called by ≥2 rounds | 84 |
| union (any round) | 104 |

A majority-of-rounds rule is count-neutral against `i = 2` (76 vs 77) but selects
a better-supported set: it drops the 7 `i = 2` calls that no majority of rounds
supports and adds the 6 above that every other round supports.

Its cost is that it rejects 17 events called at round 1 only, several at
near-whole-arm coverage (SRX10264520 1q at 1.00, SRX10831286 6p at 0.99,
SRX10831287 6p at 0.96, SRX10831280 2p at 0.94). Round agreement alone cannot
say whether those are genuine events the later rounds lost or round-1 artefacts.
Resolving them requires an orthogonal criterion — per-cell posterior support,
already computed in `results/rb_scna_probability_summary.csv` — and that check
has not yet been applied to these 17.

## 5. Recommended policy

**Superseded by §4b/§4c — retained for the reasoning, not the conclusion.** The
recommendation below assumed round 2 could serve as a fixed reference. It cannot:
round 2 is the sole dissenter in 6 cases and contributes only 3 unique calls. A
majority-of-rounds criterion is preferred, and has the further advantage of not
requiring a convergence claim at all — an event is reported when most iterations
call it, which is a stability statement rather than a fixed-point statement.

> numbat was run for four consensus iterations per sample. Because the segment
> sets do not reach a fixed point within four iterations (mean pairwise Jaccard
> 0.36 → 0.49 → 0.57, still increasing at the final step), no round can be
> designated as converged. We therefore report numbat's default second
> iteration (`i = 2`) uniformly. Round 1 was not used as the default because it
> calls 26% more segments than round 2 and its non-persisting calls are
> significantly smaller (median 6.0 vs 11.6 Mb) and lower-likelihood (median LLR
> 8.3 vs 9.9), i.e. enriched for focal artefact. Where a canonical
> retinoblastoma SCNA (1q gain, 2p gain, 6p gain, 13q loss, 16q loss) was not
> called at `i = 2`, the first iteration was inspected and the event reported if
> supported there; every such event was recovered at `i = 1`, covering a median
> 53–70% of the target arm at median LLR 75–95. Later iterations were not used
> for this purpose: iterations 3 and 4 called fewer canonical events than
> iteration 2 (76 and 71 vs 77) and recovered only 10 of the 27. Deviations are
> listed in Supplementary Table X.

**Honest caveats to keep in view.**

- 20 of 36 samples (56%) trigger the deviation. That is a majority, not a rare
  exception, and a reviewer will notice. It is defensible because the rule is
  mechanical, pre-specified, direction-fixed (always toward `i = 1`), and can
  only *add* a canonical event — but it should be stated plainly rather than
  described as occasional.
- The rule is applied only at the five canonical RB loci. Genome-wide segment
  counts and any burden statistic should come from `i = 2` alone, or round 1's
  excess focal calls will inflate them.
- This changes the current per-sample selection. The existing
  `select_numbat_round.R` picks round 1 for 19 of 36 samples on a different
  criterion and yields 98 events; a fixed `i = 2` plus the RB deviation yields
  the same 104-event union at the RB loci while keeping genome-wide counts on a
  single round.

## Files

| file | contents |
|---|---|
| `results/round_policy_report.pdf` | four-panel summary figure |
| `results/round_convergence_by_round.csv` | per sample × round: segment counts, focal counts, median size/LLR, RB events |
| `results/round1_vs_round2_segments.csv` | paired round-1 vs round-2 segment counts |
| `results/round1_segment_persistence.csv` | every round-1 segment, with whether it survives to round 2 |
| `results/round_convergence_stability.csv` | consecutive-round Jaccard, per sample per pair |
| `results/round_policy_i2_vs_union.csv` | per sample: events at `i = 2`, union across rounds, what is missed |
| `results/round_policy_missed_events.csv` | the 27 missed events with arm coverage and LLR at r2 and r1 |
| `results/round_policy_summary.csv` | the above joined to the current per-sample selection |
| `results/round_vs_r2_event_delta.csv` | per sample: events gained and lost at rounds 1, 3, 4 relative to round 2 |
| `results/rb_event_by_round_matrix.csv` | per sample x event: which rounds call it, and whether round 2 does |

Generated by `src/round_convergence_report.R`, `src/round_policy_missed_events.R`,
`src/round_convergence_stability.R`, `src/round_policy_figure.R`.
