# Hypoxia low/high splitting — strategy

*Last updated 2026-07-13. Authoritative overview of the cluster-based hypoxia
split as currently implemented. Supersedes the resolution list and marker-gate
description in [`hypoxia_split_outlier_rule.md`](hypoxia_split_outlier_rule.md),
which is retained only for the deep dive on why robust-z is used instead of a
Tukey fence.*

## Goal

Each sample's cells are partitioned into a **low-hypoxia** and a
**high-hypoxia** subset before downstream SCNA / clone analysis. The scientific
purpose is to isolate genuine hypoxic cell states so they can be characterised
separately, while keeping the low subset free of hypoxia-driven expression
structure. The split therefore has to satisfy two competing requirements:

1. **Capture real hypoxia** — clusters that are genuinely hypoxic must land in
   the high subset.
2. **Do not delete fine-grained cell-cycle-phase states** — a cluster that is
   really just a single cell-cycle phase (pure S, pure G2M, …) must stay in the
   low subset even if its hypoxia score looks elevated, because those
   phase-restricted states are exactly what the downstream analysis wants to
   keep and describe.

The tension in (2) comes from a specific confound (below), and the strategy is
built around resolving it correctly.

## Why a naive hypoxia-score threshold is not enough

Hypoxia is scored per cell (`hypoxia_score`, a module score over a hypoxia gene
set). A per-cell threshold (`subset_seu_by_expression()` on the raw score) has
two problems:

- **Cluster structure is ignored.** Hypoxia is a cell-*state*, so it should be
  called at the level of coherent groups of cells, not cell by cell.
- **The snoRNA-host-gene confound.** The hypoxia marker module includes
  snoRNA-host genes (`ZFAS1`, `GAS5`, the `SNHG` family). These are elevated in
  **proliferating S-phase cells** for reasons unrelated to hypoxia. So a pure
  S-phase (or G2M) cluster can score high on the hypoxia module and be swept
  into the "high" subset even though it is a cell-cycle state, not a hypoxic
  state. That would delete precisely the fine-grained phase clusters we want to
  keep.

So the split works on **clusters**, calls hypoxia by an **outlier rule** over
per-cluster mean scores, and then applies **two gates** that together defuse the
snoRNA-host confound.

## The strategy, in one paragraph

Sweep each sample's cells over resolutions `0.2 … 1.2`, coarse to fine. At each
resolution, flag a cluster as hypoxic when its **mean `hypoxia_score` is a
high-side statistical outlier** among that resolution's per-cluster means
(robust-z rule), keeping the flag only if the cluster **also** (a) carries at
least 2 hypoxia marker genes among its top markers (**marker gate**) and (b) is
**cell-cycle-phase-mixed** (**direct phase gate**). Every resolution is logged,
but the **drop is anchored at `r_flag`** — the *lowest* resolution at which
anything flags, i.e. the coarsest scale at which the hypoxic population first
becomes separable. All cells of every cluster flagged at `r_flag` go "high."
The survivors are then re-clustered at **`r_flag + 0.4`** and any cluster still
flagging there is dropped too (the **confirmatory clean-up**). Everything else
is "low."

This replaces an earlier rule that **unioned** flags across a fixed grid
(`0.2/0.6/1.0`) — a cell went high if it was flagged at *any* resolution. The
union was hard to reason about (the final partition was not the partition at any
single resolution) and it over-dropped.

## Step by step

### 1. Sweep resolutions; drop at the lowest STABLE one

- Resolutions `seq(0.2, 1.2, by = 0.2)`, swept coarse to fine.
- Clustering uses the `gene` assay PCA (`split_assay = "gene"`), `FindNeighbors`
  on dims 1:30 then `FindClusters` at each resolution.
- **All resolutions are logged** (`round = 1` rows in the split log) so the full
  picture stays inspectable, but only `r_flag` drives the drop. This makes the
  partition narratable: it *is* the partition at one specific resolution, plus the
  confirmatory passes.

**`r_flag` is the lowest _stable_ resolution, not the first one that flags.**
This distinction is load-bearing. A too-coarse clustering can *detect* the hypoxic
population while **under-resolving** it — capturing part of it and merging the
rest into neighbouring clusters. SRX10264518 is the worked example:

| resolution | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|
| flagged cells | **507** | 663 | 667 | 656 | 661 |

From 0.4 up, the same population reads a stable ~660 cells. Anchoring at 0.2
would leave ~24% of it in the low subset — and the confirmatory pass **cannot**
recover those cells, because once the 507-cell core is removed the remainder is
diluted and no longer forms an outlier cluster (verified: the pass flagged
nothing). So the anchor must agree with the next resolution up:

> `r_flag` = the lowest resolution whose flagged-cell count is within
> `stability_tol` (default 5%) of the next `stability_window` (default 2)
> flagging resolutions'.

For SRX10264518 this picks **0.4** (663, with 667 at 0.6 and 656 at 0.8 — all
within tolerance). For 15 of the 18 flagging samples the count is already stable
at the first flagging resolution, so the anchor is unchanged. If nothing ever
stabilises, the rule falls back to the resolution flagging the **most** cells (and
warns), so the split never knowingly under-drops.

**Why the plateau must be sustained (`stability_window = 2`).** Checking only the
*next* resolution is not enough — a monotonically rising count can contain a
spurious one-step flat spot. SRX22868102:

| resolution | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.2 |
|---|---|---|---|---|---|---|
| flagged cells | 385 | 480 | **437** | **434** | 538 | 566 |

0.6 and 0.8 agree to 0.7%, but the population is still growing — a 1-step rule
anchors at 0.6 and takes 437 cells, ~20% short of the true extent (and the capped
confirmatory pass only recovers 54 more). Requiring **two** successive agreements
rejects that flat spot (0.8 → 1.0 jumps 19%) and advances the anchor to **1.0**
(538 cells, within 4.9% of 566 at 1.2) — which matches what the union found.
- **No iterative peeling.** Peeling repeatedly re-clusters the remaining "low"
  pool and tends to chase ever-smaller marginal clusters, producing unstable,
  hard-to-interpret partitions. `n_iter` is retained in the signature only for
  call-site back-compatibility and is inert.

### 2. Outlier rule — robust z (median + k·MAD)

For each resolution, compute the per-cluster mean `hypoxia_score`. A cluster is a
score outlier when

```
mean_score > median(all cluster means) + mad_k * MAD(all cluster means)
```

with `mad_k = 3` and `MAD` scaled by 1.4826 (SD-consistent). A zero or NA MAD
(≤ 2 clusters, or many identical means) means no usable spread → nothing flagged.

**Why median + MAD and not a Tukey (Q3 + 1.5·IQR) fence:** when 2+ of a sample's
few clusters are hypoxic, those high clusters inflate Q3 and push the Tukey fence
*above* themselves, so the rule misses them. The median and MAD are dominated by
the bulk majority of non-hypoxic clusters, so the fence stays where it belongs
and the hypoxic clusters clear it. See
[`hypoxia_split_outlier_rule.md`](hypoxia_split_outlier_rule.md) for the full
argument and worked comparison.

The score outlier is **necessary but not sufficient** — it is then gated twice.

### 3. Gate 1 — marker gate (`min_hyp_markers = 2`)

The cluster must carry **≥ 2 hypoxia marker genes** among its top-20 markers
(Wilcoxon AUC, `padj < 0.5`, ranked by logFC). Marker set:

- snoRNA-host: `ZFAS1`, `GAS5`, `SNHG5/6/7/8/15/29`
- canonical hypoxia: `VEGFA`, `NDRG1`, `BNIP3`, `BNIP3L`, `SLC2A1`, `CA9`,
  `PGK1`, `LDHA`, `ENO1`, `ANKRD37`, `ADM`, `P4HA1`, `PLOD2`

This drops clusters whose high mean score is a pure artifact — e.g. a
histone/PLK2 S-phase cluster or a CENPF/MKI67 G2M cluster whose score is inflated
only by the snoRNA-host genes but which carry **no** real hypoxia markers at all.

**But the marker gate alone is not enough.** A proliferating cluster can pick up
the snoRNA-host markers (`ZFAS1`/`GAS5`/`SNHG*`) as a ride-along and thus pass the
≥2-marker test while still being a single cell-cycle phase. The marker gate
cannot tell that apart from genuine hypoxia. Hence the second gate.

### 4. Gate 2 — direct phase gate (`max_dom_frac = 0.70`)

Cell-cycle phase is scored once, up front, via `annotate_cell_cycle_without_1q()`
(Seurat `CellCycleScoring` with chr1q genes excluded so SCNA does not contaminate
the phase call), producing a per-cell `Phase`. For each cluster compute the
**dominant-phase fraction** `dom_frac = max(table(Phase)) / n_cells`.

A cluster passes the phase gate only if it is **phase-MIXED**:

```
dom_frac < 0.70          (phase_mixed = TRUE)
```

i.e. no single phase accounts for ≥ 70% of the cluster. A cluster that is ≥ 70%
one phase is **phase-restricted** and is spared from filtering even if it is a
score outlier with hypoxia markers — because that is the fine-grained phase state
we want to keep, not a hypoxic state. Clusters with no usable phase call
(`dom_frac` NA) are treated permissively (the gate does not block them; the
marker gate still governs).

This is the gate that directly enforces requirement (2): genuine hypoxia is a
metabolic state seen across cells in various phases (mixed), whereas the confound
manifests as a single-phase cluster.

**Final flag:** `is_outlier AND n_hyp_markers ≥ 2 AND phase_mixed`.

### 5. Lowest-resolution drop + confirmatory clean-up, then produce subsets

**Primary drop (`hypoxia_round = 1`).** All cells of every cluster flagged at
`r_flag` are labelled `high`.

**Confirmatory clean-up, iterated to convergence (`hypoxia_round = 2, 3, …`).**
The survivors are re-clustered at `r_flag + recluster_step` (default `+0.4`; the
sum is rounded to 1 decimal so `0.2 + 0.4` cannot become `0.6000000000000001` and
poison the column name). The same rule + both gates are applied again, and any
cluster that *still* flags is also moved to `high`. The pass then **repeats** at
`r_flag + 2·step`, `r_flag + 3·step`, … **until a pass flags nothing.**

A single confirmatory pass is **not** sufficient. Empirically, SRX10264518 flags
at `r_flag = 0.2`, is clean at 0.6 — and yet the old union rule caught a further
~160-cell hypoxic pocket that only separates at resolution 1.0. Stopping after
one pass would leave that residual in the low subset, which is the input to the
SCNA/clone analysis. Iterating to quiescence closes that hole: the loop only
stops once the retained cells contain no flaggable hypoxia at any resolution it
has reached.

Each pass gets its own `hypoxia_round` index, so the log shows **which
resolution finally separated each pocket**.

**The ladder is capped at `max(resolutions)`** (1.2). It must not climb past the
swept range: above it, clusters fragment and the MAD fence becomes easier to
clear, so an unbounded ladder *manufactures* drops at resolutions nobody asked to
search. SRX10264520 is the cautionary case — it flags nothing through 1.0, flags
74 cells at 1.2, and an uncapped ladder then took a further 105 cells at 1.6, on a
sample that the production rule left entirely alone. The cap is why
`max_recluster_res` defaults to `max(resolutions)` rather than a free constant.
The trade-off is explicit: a sample whose `r_flag` is already 1.2 has no room for
a confirmatory pass and gets the primary drop only.

The iteration also terminates if the pool drops below 20 cells, or after
`max_confirmatory_rounds` (default 5, which **warns** — that means it did not
converge and residual hypoxia may remain). Set `confirmatory_drop = FALSE` to run
a single pass that records the residual in the log without dropping it.

Everything not dropped is `low` (`hypoxia_partition` metadata column). The
labelled object is written to `*_seu_hypoxia_labeled.rds`, then the established
`subset_seu_by_expression()` path re-clusters, re-scores markers, recomputes UMAP
and hash/barcode metadata for each subset, emitting:

- `*_hypoxia_low_seu.rds` (low subset, `gene` assay)
- `*_hypoxia_high_seu.rds` (high subset, `SCT` assay)

so every downstream target consumes the same file shapes as before. Note the
`r_flag + 0.4` clustering is **diagnostic only** — the persisted low object
carries `subset_seu_by_expression()`'s own full resolution grid
(`seq(0.2, 1, 0.2)`), which is what downstream targets key on
(`SCT_snn_res.0.6`).

### Known cost — honest framing

Anchoring at the lowest flagging resolution is **not** biologically superior to
the union on requirement (2). `r_flag` is by construction the *least-separated*
scale, so a cluster flagged there is the coarsest, most-merged version of the
hypoxic population — dropping it is if anything more likely, not less, to take
neighbouring cells with it. On "don't delete fine-grained phase states" the new
rule is **neutral** to the union, not better.

The genuine risk of anchoring low was **under-dropping**: residual hypoxia that
the union would have caught at a higher resolution surviving into the low subset.
That is the expensive error direction, and it is not hypothetical — it showed up
on the first sample tested (SRX10264518, above). The **iterated** confirmatory
pass is what neutralises it: the loop keeps re-clustering the survivors at finer
and finer resolution until nothing flags, so a pocket that only separates at 1.0
is still caught, just as it would have been under the union.

What remains, then, is a rule that drops **no more** than the union (the primary
drop is anchored at the coarsest separable scale) and — once converged — leaves
**no more residual** than the union either. The wins are **interpretability** (the
partition is a real clustering plus a stated sequence of drops, not a set union),
a **more surgical** drop, and a log that records **which resolution separated each
pocket** instead of collapsing everything into one undifferentiated union.

## What the two gates buy — cohort evidence

Offline diagnostic (`src/eval_hypoxia_cluster_phase.R`, resolutions 0.2/0.6/1.0)
scored all **92 ungated score-outlier clusters** across 25 samples and classified
each by phase purity + hypoxia-marker content:

| Diagnostic call                     | n  | Two-gate decision |
|-------------------------------------|----|-------------------|
| genuine hypoxia, phase-mixed        | 28 | **KEPT (filtered to high)** |
| hypoxia markers but phase-restricted| 32 | spared → stays low |
| phase artifact (no real markers)    | 19 | spared → stays low |
| ambiguous / non-hypoxia             | 13 | spared → stays low |

The two gates **keep exactly the 28 genuine, phase-mixed hypoxic clusters** and
**spare the other 64** — critically including the 32 phase-restricted clusters
that a **marker-only** rule would have wrongly filtered. 4 of the 16 res-0.2
clusters also show phase-splitting at higher resolution, the failure mode the
strategy is designed to avoid.

The production split (with the gates on) flagged **31 clusters across 18 samples**
in `hypoxia_split_log_all.csv`; **every** flagged cluster is `phase_mixed = TRUE`,
and the maximum `dom_frac` among them is **0.696** (< 0.70), confirming the phase
gate is binding. (31 vs 28 is Louvain seed jitter between the diagnostic and
production runs, not a rule difference.)

## Outputs and where to look

- **Per-cell labels:** `hypoxia_partition` (`low`/`high`) on the Seurat object;
  drives `*_hypoxia_low_seu.rds` / `*_hypoxia_high_seu.rds`.
- **Per-sample decision log:** `results/hypoxia_cluster_split/<SRX>_hypoxia_split_log.csv`
  — one row per (round × resolution × cluster) with `n_cells`, `top_markers`,
  `n_hyp_markers`, `matched_genes`, `mean_score`, `dom_frac`, `outlier_fence`,
  `is_outlier`, `phase_mixed`, `flagged`. `round = 1` rows are the full
  coarse-to-fine sweep (every resolution, whether or not it drove the drop);
  `round = 2, 3, …` rows are the successive confirmatory passes on the survivors
  at `r_flag + k·0.4`. The last confirmatory round present is the one that
  converged (flagged nothing).
- **Collated log:** `results/hypoxia_cluster_split/hypoxia_split_log_all.csv`
  (`hypoxia_split_log_collated` target) — all samples stacked; rclones to
  `gdrive:.../hypoxia_rebuilt/hypoxia_split_log/`.
- **Per-sample QC PDF:** `<SRX>_hypoxia_split.pdf` — UMAP by partition, UMAP by
  round, partition marker dotplot.
- **Diagnostic stage collages** (`write_diagnostics = TRUE`, only for samples that
  flag), in `results/`, each a marker-heatmap / phase-scatter patchwork grouped on
  the split clustering:
  1. `<SRX>_hypoxia_rflag<R>_before_seu.rds__filtered_heatmap_phase_scatter_patchwork.pdf`
     — at `r_flag`, **before** the drop (the flagged clusters are visible).
  2. `<SRX>_hypoxia_rflag<R>_after_seu.rds__filtered_heatmap_phase_scatter_patchwork.pdf`
     — same clustering, **after** the drop.
  3. `<SRX>_hypoxia_recluster<R'>_seu.rds__filtered_heatmap_phase_scatter_patchwork.pdf`
     — **one per confirmatory pass**, at each `R' = r_flag + k·0.4`, showing that
     pass's clustering with any flagged cluster still visible *before* it is
     removed. The final one is the converged state (nothing flagged), so reading
     them in order shows the residual being peeled away until the data is clean.

## Where it lives in code

- `numbat_helpers/R/hypoxia_cluster_split.R`
  - `identify_hypoxia_clusters()` — the per-resolution rule + gates + decision table.
  - `split_hypoxia_by_clusters()` — the sweep, `r_flag` drop, confirmatory pass,
    phase scoring, subset emission.
  - `.write_hypoxia_stage_heatmap()` — the three diagnostic stage collages.
- `R/pipeline_targets_seurat.R`
  - `hypoxia_partition_paths` — calls `split_hypoxia_by_clusters(..., max_dom_frac = 0.70)`.
    Note this target has `cue = tar_cue(command = FALSE)`, so editing its command
    does **not** trigger a rebuild — `tar_invalidate(hypoxia_partition_paths)` first.
  - `hypoxia_split_log_collated` — stacks the per-sample logs.
- `R/pipeline_targets_figures.R`
  - `hypoxia_rebuilt_gdrive` — rclones the low-hypoxia deliverables + the split log.

## Key parameters (defaults)

| Parameter          | Value              | Meaning |
|--------------------|--------------------|---------|
| `resolutions`      | `seq(0.2, 1.2, 0.2)` | Louvain resolutions swept coarse→fine; all logged, the lowest **stable** one (`r_flag`) drives the drop |
| `stability_tol`    | `0.05`             | `r_flag` = lowest resolution whose flagged-cell count is within 5% of the next resolutions' |
| `stability_window` | `2`                | how many successive resolutions the anchor must agree with (1 admits spurious flat spots) |
| `recluster_step`   | `0.4`              | confirmatory pass `k` re-clusters survivors at `r_flag + k·0.4` |
| `confirmatory_drop`| `TRUE`             | drop clusters still flagging in a confirmatory pass, iterating until one flags nothing |
| `max_confirmatory_rounds` | `5`         | safety cap on confirmatory passes; hitting it **warns** (did not converge) |
| `max_recluster_res`| `max(resolutions)` = `1.2` | ladder ceiling — never search past the swept range |
| `write_diagnostics`| `TRUE`             | emit the stage collages |
| `n_iter`           | `1`                | inert; kept for call-site back-compatibility |
| `mad_k`            | `3`                | robust-z multiplier for the outlier fence |
| `min_hyp_markers`  | `2`                | marker gate: ≥ 2 hypoxia markers in top-20 |
| `max_dom_frac`     | `0.70`             | phase gate: flag only if dominant-phase fraction < 0.70 |
| `top_n`            | `20`               | marker window scanned for the hypoxia overlap |
| `split_assay`      | `"gene"`           | assay for partition clustering |
| `low_assay`/`high_assay` | `"gene"`/`"SCT"` | subset re-clustering assays |

Set `min_hyp_markers = 0, max_dom_frac = NULL` to recover the original
score-only rule (used by the offline phase diagnostic).

## Related docs

- [`hypoxia_split_outlier_rule.md`](hypoxia_split_outlier_rule.md) — deep dive on
  the robust-z vs Tukey choice.
- `src/eval_hypoxia_cluster_phase.R` — the offline per-cluster phase diagnostic.
