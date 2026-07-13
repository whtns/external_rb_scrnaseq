# Hypoxia low/high split — outlier mean-score rule

> **Partially superseded.** For the current end-to-end strategy — including the
> two gates (marker + direct phase) and resolutions `0.2/0.6/1.0` — see
> [`hypoxia_splitting_strategy.md`](hypoxia_splitting_strategy.md). This doc is
> retained for its deep dive on the robust-z (median + k·MAD) vs Tukey-fence
> choice; its resolution list and "marker gate removed" paragraph below are out
> of date.

## Summary

The hypoxia cluster split partitions each sample's cells into a **low-hypoxia**
and a **high-hypoxia** subset before downstream SCNA analysis. As of this
revision the rule is:

- **One round**, no iterative peeling (`n_iter = 1`).
- Cells are clustered at **multiple resolutions** (`0.2`, `0.4`, `0.6`). A whole
  cluster is moved to the high subset when its **mean `hypoxia_score` is a
  high-side statistical outlier** among that resolution's per-cluster means,
  using the **robust-z (median + 3·MAD)** rule below.
- Hits are **unioned across resolutions**: a cell is "high" if it lands in a
  flagged cluster at any of the three resolutions.
- Flagging is **independent of hypoxia marker-gene presence**. The previous
  scheme also flagged a cluster when ≥2 of its top markers were hypoxia genes
  (`marker_hit`); that rule has been removed. Marker information is still
  reported in the log for interpretation but does not drive the split.

Implementation: `numbatHelpers::split_hypoxia_by_clusters()` and
`numbatHelpers::identify_hypoxia_clusters()` in
`numbat_helpers/R/hypoxia_cluster_split.R`. Wired via the
`hypoxia_partition_paths` target in `R/pipeline_targets_seurat.R`.

## Outlier rule (chosen): robust z / median + k·MAD

For each resolution, compute the per-cluster mean `hypoxia_score`. A cluster is
flagged when

```
mean_score > median(means) + 3 * mad(means)
```

`mad()` is the median absolute deviation, scaled by 1.4826 so it is consistent
with the SD for normal data; `k = 3` ≈ a 3-SD cutoff. Because the median and MAD
are dominated by the **bulk majority** of (non-hypoxic) clusters, the threshold
is **not** pulled up by the elevated clusters themselves — so it still catches
**multiple** hypoxic clusters in the same sample. Degenerate cases (`MAD = 0`,
e.g. ≤2 clusters or many identical means) flag nothing, which is the desired
conservative behavior — a sample with no outlier cluster simply gets an empty
high subset.

The threshold is recorded per row as `outlier_fence`, and the decision as
`is_outlier` (equal to `flagged`).

### Why not the Tukey fence

The classic boxplot rule `mean_score > Q3 + 1.5·IQR` was tried first and
**failed validation**: on the reference samples, when 2 of a sample's ~5
clusters (or the single hypoxic cluster among only 3) are hypoxic, those high
clusters sit in the top quartile and **inflate Q3**, pushing the fence above
themselves. It correctly flagged only 1 of 4 reference samples
(SRX10264517/18/20 missed; SRX10264519 caught — the one with a single elevated
cluster). The MAD rule reproduced all four expected labelings
(SRX10264517→{2,4}, SRX10264518→{4,5}, SRX10264519→{2}, SRX10264520→{2}).

## Alternative rules (not used)

These can be swapped in by editing the threshold computation in
`identify_hypoxia_clusters()`:

1. **Tukey upper fence**

   ```
   mean_score > Q3 + 1.5 * IQR        (IQR = Q3 - Q1, over the per-cluster means)
   ```

   Classic non-parametric boxplot outlier. Works when a single outlier sits
   among many tightly-clustered values, but Q3 is inflated when several clusters
   are elevated (see "Why not the Tukey fence" above).

2. **Mean + k·SD**

   ```
   mean_score > mean(means) + 2 * sd(means)
   ```

   Parametric Gaussian rule. Simple, but the **outlier itself inflates the mean
   and SD**, which can mask the very cluster it is meant to catch — especially
   with a small number of clusters. `k = 2` flags roughly the top 2.5% under
   normality.

## Per-sample log columns

Written to `results/hypoxia_cluster_split/<SRX>_hypoxia_split_log.csv` and
collated into `results/hypoxia_cluster_split/hypoxia_split_log_all.csv`
(one row per `round` × `resolution` × `cluster`):

| column          | meaning                                                        |
|-----------------|----------------------------------------------------------------|
| `sample_id`     | sample (SRX/SRR id)                                            |
| `round`         | peeling round (always `1` in the single-round scheme)         |
| `resolution`    | Louvain resolution for this clustering (`0.2`/`0.4`/`0.6`)     |
| `pool_n`        | cells fed to this round's clustering                          |
| `cluster`       | cluster id at this resolution                                 |
| `n_cells`       | cells in this cluster                                         |
| `top_markers`   | top 5 marker genes (by logFC, padj < 0.5) — interpretation     |
| `n_hyp_markers` | how many of the top `top_n` markers are hypoxia genes (info)   |
| `matched_genes` | which hypoxia genes matched (info)                            |
| `mean_score`    | mean `hypoxia_score` of the cluster                           |
| `outlier_fence` | Tukey upper fence for this resolution's per-cluster means      |
| `is_outlier`    | `mean_score > outlier_fence`                                  |
| `flagged`       | cluster moved to the high subset (= `is_outlier`)             |
