# Hypoxia vs mitochondrial score boxplots

**Target:** `hypoxia_mito_score_boxplots`
**Output:** `results/hypoxia_cluster_split/hypoxia_mito_score_boxplots.pdf` (2 pages)
**Drive:** `gdrive:rb_scrnaseq/hypoxia_rebuilt/mito_score_boxplot/{paper_retained,other}/`
**Build:** `sbatch build_hypoxia_mito_boxplots.sbatch`

## Why

The hypoxia split flags a cluster when its mean `hypoxia_score` clears a
`median + 3*MAD` fence. That score is a **composite**, built by
`add_hypoxia_score()`:

```r
seu$MT           <- seu$MT * -1                       # NEGATED
seu$hypoxia_score <- rescale(rowMeans(hypoxia, MT), c(0, 1))
```

So a cluster can clear the fence two different ways:

1. genuinely high HALLMARK_HYPOXIA expression, or
2. merely **low mitochondrial content** — which the negation turns into a high
   score just as effectively.

The existing boxplot
(`src/plot_hypoxia_mean_score_boxplots.R` → `hypoxia_mean_score_boxplots.pdf`)
plots the composite, which is the right thing to show for the *decision* but
cannot distinguish those two cases. This figure shows the two ingredients
separately, so a reader can tell which half drove each flag.

On the smoke samples the answer was already visible: flagged and spared clusters
sit **above** their sample's bulk on the hypoxia page and **far below** it on the
mitochondrial page (≈0.5 against a bulk of ~1.5 for SRX10264517/18) — the
mitochondrial half is doing a large share of the flagging work.

## Finding: `hypoxia_score` is ~99.8% an inverted mitochondrial score

`src/diag_hypoxia_mito_correlation.R` (33 samples, per cell,
`results/hypoxia_cluster_split/hypoxia_mito_correlation.csv`):

| quantity | median | range |
|---|---|---|
| `r(hypoxia_score, mito_score)` | **-0.998** | -0.975 … -0.999 |
| `r(hypoxia_score, hypoxia)` | **0.387** | 0.067 … 0.572 |
| `r(hypoxia, mito_score)` raw modules | -0.312 | -0.562 … 0.157 (32/33 negative) |
| mito share of ingredient sd | **0.935** | 0.818 … 0.962 |

The composite is `mean(hypoxia, -mito)`, an **unweighted** mean of two scores
whose spreads differ by roughly 14x: the MT module is 5 very highly expressed
genes, HALLMARK_HYPOXIA is ~200 moderately expressed ones. The mitochondrial
half therefore carries ~93% of the ingredient standard deviation, and the mean
collapses to `-mito/2`.

Consequences:

- The statistic the split thresholds on is, to three decimal places, the
  negative of the mitochondrial score. "Flagged as hypoxic" reads in practice as
  **"low mitochondrial content"**.
- The hypoxia gene module contributes weakly (r = 0.39 with the composite).
- The raw modules *are* genuinely anti-correlated (median r = -0.31, 32/33
  samples), so some of the inversion is real biology — but nothing like enough
  to explain -0.998. The rest is the unweighted average.
- This is consistent with what the boxplot shows: flagged clusters sit modestly
  above their sample's bulk on the hypoxia page and far below it on the
  mitochondrial page.

The obvious remedy is to standardise each ingredient (z-score within sample)
before averaging, so the two halves contribute equally. **Not done here** —
changing the score changes every split decision and invalidates the whole
hypoxia chain, so it is the user's call.

## Layout

Two pages, one per score, each with the familiar layout: one facet per sample,
clustering resolution on x, one point per cluster, coloured by the three-state
split decision (`not an outlier` / `outlier, spared by gates` / `flagged → high
subset`) and outliers labelled `c<id> (n=<cells>)`.

Two deliberate departures from the composite plot:

- **One page per score, not two dodged boxes in one panel.** The mitochondrial
  module is 5 very highly expressed genes and lands around 0.5–2.5;
  HALLMARK_HYPOXIA is ~200 genes and lands near 0.05. Dodged into one panel on a
  shared axis, the hypoxia boxes flatten to a line. (This was tried first and
  rejected on inspection.)
- **Free y per sample facet.** The composite plot can share one axis across
  samples because `hypoxia_score` is rescaled to [0,1] *within each sample*.
  These are raw module scores, whose absolute level tracks library depth and
  composition, so a shared axis is neither comparable nor readable — one
  high-range sample flattens the other 32. **Compare within a panel, not across
  panels.**

`mean_mito` is stored un-negated, so **higher = more mitochondrial** — the
opposite orientation to the `MT` column on the Seurat objects.

## Where the numbers come from

Two paths, preferred first.

### 1. The split log (preferred, forward-looking)

`identify_hypoxia_clusters()` now records `mean_hypoxia` and `mean_mito`
alongside `mean_score` in its per-cluster detail frame, so every
`split_hypoxia_by_clusters()` run from that change onward writes them into
`results/hypoxia_cluster_split/<sample>_hypoxia_split_log.csv` and thence into
the collated `hypoxia_split_log_all.csv`. They are recorded for interpretation
only — no gate or fence reads them.

### 2. Replay (fallback, used today)

The existing log predates those columns, and `hypoxia_partition_paths` is pinned
with `cue(command = FALSE)`, so editing the split does **not** re-run it. Until
someone deliberately re-runs the split, `replay_cluster_score_means()` recovers
the two means by reproducing the split's round-1 sweep from the
`*_seu_hypoxia.rds` objects: same assay, same
`FindNeighbors(dims = 1:30, k.param = min(20, n - 1))` on the stored `pca`, same
`FindClusters` per resolution. Only `presto::wilcoxauc()` is skipped — it feeds
the marker gate, not the score means.

**The replay is verified, not trusted.** Every replayed cluster is joined to the
logged one on `(sample_id, resolution, cluster)` and kept only if `n_cells` is
identical *and* `mean_score` matches within `1e-6`; together those pin the
cluster assignment. If nothing reconciles the function refuses to plot rather
than showing means that might belong to a different clustering.

Measured on the 3-sample smoke test: **191/191 clusters reconciled, max
`mean_score` difference 5.55e-17** — i.e. the replay reproduces the split's
clustering exactly, to floating-point noise.

Only **round 1** is plotted. The confirmatory rounds re-cluster a shrinking pool
at resolutions computed as `r_flag + k*recluster_step`, so their cluster ids do
not share a grid with the sweep; restricting to round 1 also makes the log and
replay paths cover exactly the same rows.

## Files

| File | Role |
|---|---|
| `numbat_helpers/R/hypoxia_mito_score_boxplots.R` | `replay_cluster_score_means()`, `plot_hypoxia_mito_score_boxplots()` |
| `numbat_helpers/R/hypoxia_cluster_split.R` | `identify_hypoxia_clusters()` now logs `mean_hypoxia` / `mean_mito` |
| `R/pipeline_targets_figures.R` | `hypoxia_mito_score_boxplots` target; `mito_score_boxplot` entry in `hypoxia_rebuilt_gdrive` |
| `src/smoke_hypoxia_mito_boxplots.R` + `smoke_hypoxia_mito_boxplots.sbatch` | 3-sample smoke test; prints the reconciliation table |
| `build_hypoxia_mito_boxplots.sbatch` | production build + Drive sync |

## Notes

- The build script invalidates the target explicitly: `targets` does not track
  package function bodies, so editing the plotting function does not invalidate
  it on its own.
- The base `pdf()` device emits `mbcsToSbcs` warnings substituting `-` for `—`
  and `->` for `→`. Cosmetic, and the same happens in the composite-score plot.
