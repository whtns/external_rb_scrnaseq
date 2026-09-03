# Sequencing depth: the panel, and why there is no new cutoff

**Issue:** [#43](https://github.com/whtns/external_rb_scrnaseq/issues/43)
**Panel:** `numbatHelpers::plot_cell_depth_panel()`, one row in each
`results/<sample>_summary.pdf`
**Diagnostic:** `src/diag_cell_depth_provenance.R` → `results/cell_depth_provenance.csv`
**Run:** `sbatch diag_cell_depth.sbatch`

Issue #43 asked to "plot read count on sample summaries and possibly exclude low
read count cells". The panel is delivered. The exclusion is **not**, and this
document is the evidence for why.

## What the data says

From `batch_hashes.sqlite` / `cell_qc_values` — 261,618 cells across 38 samples,
at the unfiltered stage.

### 1. A read-count filter already exists

`filter_sample_qc()` (`numbat_helpers/R/plot_functions_3.R:404`) and the
annotation filter in `generate_filtering_cell_counts()`
(`numbat_helpers/R/targets_helpers.R:25`) both apply:

```r
percent.mt < 10 & nCount_gene > 1000 & nFeature_gene > 1000
```

### 2. But numbat never sees it

`pipeline/scripts/run_numbat.R:119-147` takes its cell set from the **unfiltered**
`*_seu.rds` and subsets on cell **type** only. numbat's own `min_depth` defaults
to **0** and is not overridden.

So CNV inference runs on cells the Seurat analysis later discards. That is the
single most actionable finding here, and it is what the panel is drawn to show —
the caption says *"numbat is run on these cells, before that filter."*

### 3. A cohort-wide cutoff is not supported

| threshold | cells below | share |
|---|---|---|
| 500 | 0 | 0.00% |
| 1,000 | 8,125 | 3.11% |
| 2,000 | 25,392 | 9.71% |

Nothing sits below 500, and only 3% below the filter the pipeline already
applies. There is no population of junk cells to remove cohort-wide.

### 4. The existing floor is *not uniform* — this is the real problem

Sample minima fall into two clear regimes:

| regime | samples |
|---|---|
| floor ≈ 500 | 16 |
| floor ≈ 1250–1530 | 22 |

and the split tracks when the object was processed, not the biology:

| regime | `_seu.rds` mtime | samples |
|---|---|---|
| high (≥1250) | 2024-09-26 | 1 |
| high (≥1250) | 2026-06-01 | 6 |
| high (≥1250) | 2026-06-02 | 9 |
| high (≥1250) | 2026-06-03 | 6 |
| low (~500) | 2024-09-26 | 9 |
| low (~500) | 2026-03-24 | 6 |
| low (~500) | 2026-06-03 | 1 |

The 2026-06-01/02 batch is uniformly high-floor; 2026-03-24 is uniformly
low-floor. This is a **processing-vintage confound**, and it is what actually
wants fixing — not a new threshold layered on top of an inconsistent one.

### 5. Depth does not predict numbat outcome — except at the extreme

Median depth spans **22.8-fold** across samples (831 → 18,980), yet
`spearman(median depth, n_rb_events) = 0.282` — weak.

Exactly one sample has incomplete numbat artefacts, and it is the lowest-depth
sample in the cohort:

| sample | median | % <1000 | depth rank | RB events | missing artefacts |
|---|---|---|---|---|---|
| SRX10031191 | 831 | 59.8% | 1 / 35 | 0 | `clone_post;gtree;mut_graph;treeML` |

The next-lowest sample has median 2,501 and 0% below 1,000. So the depth/quality
association is carried entirely by this one sample. Its panel shows a clearly
**bimodal** distribution: a large mode of ~870 cells below the filter line and a
small proper mode near 10–20k. numbat ingested all of it.

## Conclusion

- **Delivered:** the depth panel, on every sample summary, with the pipeline's own
  `nCount_gene > 1000` line marked.
- **Not delivered, deliberately:** a new low-count exclusion. The data does not
  support one, and adding a second threshold on top of an inconsistent first one
  would make provenance harder to reason about, not easier.
- **Recommended next**, as separate decisions:
  1. Harmonise the upstream floor so the cohort is filtered once, consistently.
  2. Decide whether numbat should run on the filtered cell set rather than the
     unfiltered one — currently it does not, and nothing says so anywhere.
  3. Treat SRX10031191 as a known-bad sample on depth grounds rather than
     hoping a cohort threshold catches it.

Each of those changes CNV calls, so none was made here.

## Implementation notes

The panel reads `cell_qc_values` from `batch_hashes.sqlite` via
`get_sample_cell_depth()`, using `connect_hash_db()` / `db_retry()` per
AGENTS.md. Re-reading the Seurat objects for the same numbers would cost minutes
per sample for no gain.

It is added as **one row beneath the three columns**, not one per column: it
describes the unfiltered cell set, which is shared, so repeating it three times
would say the same thing three times.
