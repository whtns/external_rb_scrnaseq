# Does genome dosage explain the depth–clone association?

**Resolves the direction question left open by [depth_clone_association.md](depth_clone_association.md).**

- Script: `src/diag_clone_dosage_depth.R` · Job: `diag_clone_dosage_depth.sbatch` (11784845, 4m42s, 1.6 GB)
- Outputs: `results/clone_dosage_depth_samples.csv`, `results/clone_dosage_depth_clones.csv`,
  `results/clone_dosage_phi_check.csv`, `results/clone_dosage_depth.pdf`
- 36 samples, 198 clones. Read-only: reads `output/numbat_sridhar/*_numbat.rds` and
  `cell_qc_values`; writes only `results/`. No numbat rerun, no `_targets_r431` access.

## The question

Part 1 found that tumor clones separate on cell depth in 30/30 testable samples
(median epsilon-squared 0.047), and I could not say which way the arrow pointed.
The premise I had used to propose that test was wrong: I claimed sibling
subclones "have no biological reason to differ in depth," but subclones are
*defined* by gains and losses, and genome dosage genuinely changes RNA content.
A gain-bearing subclone should have more UMIs, and that would produce exactly
the observed signal without any artifact.

Dosage is testable because it makes a **quantitative** prediction, not just a
sign. Each clone's `GT_opt` names its `seg_cons` segments; `segs_consensus`
gives each a `cnv_state`, a length, and (via `gtf`) a gene content. So:

```
predicted log2 depth fold = log2(1 + Σ w_s · δ_s)
  δ : amp +0.5, bamp +1.0, del −0.5, bdel −1.0, loh 0   (copies/2 − 1)
  w : that segment's share of the assessed genome, by gene count (primary)
      and by bp length (secondary)
```

`GT_opt` tokens matched `seg_cons` at **100% in every sample**, so nothing is
being silently scored as zero.

## 0. The dosage model checks out against numbat's own measurements

`segs_consensus$phi_mle` is numbat's measured expression fold per segment. It is
not used to build the prediction, so it is an independent check on the δ values:

| state | assumed | measured median phi | IQR | n |
|---|---|---|---|---|
| amp | 1.5 | **1.49** | 1.31–1.71 | 232 |
| del | 0.5 | **0.62** | 0.55–0.74 | 368 |
| loh | 1.0 | **1.00** | 0.92–1.09 | 141 |
| bamp | 2.0 | 1.63 | 1.48–1.72 | 10 |
| bdel | 0.0 | 0.59 | 0.53–0.62 | 20 |

The two states that carry the analysis — `amp` and `del`, 600 of 771 segments —
land essentially on the assumed values. `bamp` and `bdel` are overestimated by
the model, but they are 30 segments total, and **the error runs in the direction
that favours dosage**: substituting measured phi would shrink predicted effects
and widen the gap reported below. The prediction is therefore conservative.

## 1. Dosage is real — gains genuinely make clones deeper

| test | median rho | positive | sign test |
|---|---|---|---|
| clone-level, all clones | +0.520 | 22/30 | p = 0.016 |
| clone-level, tumor only | +0.555 | 20/25 | p = 0.0041 |
| cell-level, all cells | +0.088 | 26/35 | p = 0.0060 |
| **cell-level, tumor only** | **+0.104** | **24/30** | **p = 0.0014** |
| cell-level, tumor, bp-weighted | +0.092 | 23/30 | p = 0.0052 |

So the biological explanation is not a fiction. Clones carrying more genome are
measurably deeper, consistently across the cohort, and the result is stable to
swapping the weighting scheme.

## 2. But dosage is far too small to account for the observed spread

| comparison | observed | predicted | ratio |
|---|---|---|---|
| normal-referenced max abs log2 fold | 1.268 | 0.024 | **35×** |
| **tumor-only spread across clones** | **0.723** | **0.026** | **41×** |

`lm(observed ~ predicted)` over 162 tumor clones: **slope +10.9** (SE 2.38,
p = 1e−05), **R² = 0.115**. A slope of 1 would mean dosage accounts for the
depth differences.

The tumor-only row is the one that matters — it excludes the normal clone
entirely, so the large normal-vs-tumor depth confound (part 1: median log2 ratio
+0.75) cannot contribute. **Observed spread exceeds the dosage prediction in
30/30 samples.**

The reason is arithmetic: gains and losses largely cancel in net dosage. Within-
sample dosage range has a median of just **0.020** — a predicted depth
difference of **2%** between the most and least dosed clone. Observed depth
differences between those same clones are tens of percent, and in the extreme
cases (SRX10031193, SRX10031194, SRX14116947) 4–8×.

## What this establishes

**Dosage is a genuine but minor contributor.** It explains ~11% of the variance
in clone depth differences and roughly 1/41 of their magnitude. The depth–clone
association reported in part 1 is therefore **not** mainly a dosage artifact of
my own test design — the alternative explanation I raised against myself turns
out to be real but small.

**What it does not establish.** "Not dosage" is not the same as "artifact."
Subclones can differ in depth for other genuine biological reasons that this
analysis does not model — proliferative state, cell size, differentiation stage,
total RNA content independent of copy number. Ruling out the specific dosage
mechanism narrows the field; it does not single out read depth as the cause.

**The causal test is still downsampling**: take the deep samples, downsample
reads to the shallow ones' depth, rerun numbat, and see whether the clone
partition survives. That is the only design that breaks the correlation by
intervention. It requires numbat reruns and so must respect the
`SRX11133592/93/94` hold-out.

## Practical consequence, unchanged from part 1

The solid, immediately actionable finding remains the part 1 one:
`spearman(nCount_gene, p_opt) > 0` in **36/36** samples, median **+0.319**.
`p_opt` is a posterior confidence, so unlike depth-vs-clone it has no competing
biological reading. Cells in the shallow tail carry clone labels numbat itself
is unsure of, and nothing downstream propagates that uncertainty.
