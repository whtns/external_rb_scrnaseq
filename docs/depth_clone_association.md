# Are numbat clones an artifact of read depth?

Per-cell analysis, 36 samples / 225,525 cells. Job 11750635.
Script `src/diag_depth_clone_association.R`; outputs
`results/depth_clone_association{.csv,.pdf}`, `depth_clone_per_clone.csv`,
`depth_clone_cells.csv.gz`.

## Why this exists

Everything #43 delivered is between-sample and descriptive: a marginal
`nCount_gene` panel per summary, plus one weak sample-level correlation
(`spearman(median depth, n_rb_events) = 0.282`, p = 0.10, n = 35; falling to
0.219 / p = 0.21 once the single broken sample SRX10031191 is dropped). None of
it can address whether clones are depth artifacts, because that is a
within-sample, per-cell question.

`nb$clone_post` and `cell_qc_values` share barcodes exactly (4085/4085 on the
pilot; 34/36 samples join at 100%, the two exceptions being SRX11133592 and
SRX11133594 where numbat carries 10,000 cells against ~9,000 in the QC table).

## Design, and the control that matters

numbat's clone 1 has an empty `GT_opt` — it is the normal/diploid clone. So a
depth split across *all* clones is expected and uninformative. The reported test
is therefore computed twice: over all clones, and over **tumor clones only**.
Effect size is epsilon-squared from Kruskal-Wallis, because at n in the
thousands every p-value is significant and only the effect size is readable.

## Results

### 1. The normal/tumor depth confound is large

Tumor cells are deeper than normal cells in **33 / 36** samples, median
`log2(tumor/normal) = +0.75` — tumor cells carry ~1.7x the counts.

This is the most consequential single number here. The normal-vs-malignant call
is what most downstream analysis leans on, and it is exactly the axis most
entangled with depth.

### 2. Tumor subclones separate on depth in every sample, but weakly

| | value |
|---|---|
| samples testable (>= 2 tumor clones, >= 50 cells) | 30 |
| BH q < 0.05 | **30 / 30** |
| median epsilon-squared | **0.047** |
| samples with epsilon-squared > 0.06 ("moderate") | 11 |
| range | 0.0014 - 0.249 |

Universal significance, mostly small effects: clone explains ~5% of depth
variance in the median sample. Not an artifact of clone count
(`rho(n_tumor_clones, eps2) = +0.35, p = 0.06`) or of cell count
(`rho = +0.08, p = 0.67`). It holds in the five samples with only two tumor
clones, where there is least room to fish.

### 3. numbat is less certain about shallow cells — in all 36 samples

`spearman(nCount_gene, p_opt)`: **positive in 36 / 36 samples**, median
**+0.319**, quartiles 0.20-0.46, max 0.73.

`spearman(nCount_gene, mean p_cnv)` over called segments: positive in 25/36,
median +0.144.

Neither tracks sample median depth (`rho = -0.21, p = 0.22`), so this is not a
property of shallow *samples* — it is a property of shallow *cells*, everywhere.

## What this does and does not establish

**Established.** Clone assignment is systematically less confident for
low-depth cells, in every sample without exception, and often strongly.
`p_opt` is a posterior confidence, not a dosage quantity, so there is no
biological reading of this: it is loss of inferential power at low counts.
Cells in the shallow tail have clone labels that numbat itself is unsure of, and
nothing downstream currently propagates that uncertainty.

**Not established — and analysis 1 cannot establish it.** Whether the depth/clone
*association* means depth drives the partition. The premise this test was
designed on ("sibling subclones have no biological reason to differ in depth")
is weaker than it looked. Subclones are defined by copy-number gains and losses,
and genome dosage genuinely changes RNA content: a gain-bearing subclone should
have more counts. That mechanism produces exactly the observed signal and is not
an artifact. The direction is unresolved.

## Next step that would resolve the direction

Per tumor clone, map the segments in `GT_opt` to their `cnv_state` in
`segs_consensus`, compute net genome dosage, and correlate that with the clone's
median depth within each sample. If depth tracks dosage, the association in
section 2 is biology. If it does not, depth is driving the partition. Cheap:
same objects, no numbat rerun.

Separately, the only causal test of the detection question remains downsampling
a few high-depth samples and rerunning numbat, where the truth is fixed by
construction.

## Caveats

- 36 of 39 samples: SRX10031191 has no `clone_post` (it is the broken low-depth
  sample from #43), and two more lack usable joins.
- `nCount_gene` comes from `cell_qc_values`, i.e. the Seurat-side count. numbat
  ingests the unfiltered cell set with `min_depth = 0` (see
  `docs/cell_depth_provenance.md`), so the shallow cells analysed here are cells
  numbat actually used.
- Kruskal-Wallis tests location shift only; clones differing in depth *spread*
  but not median would not register.

---

## Follow-up: the direction question is now resolved

The "directionally unresolved" caveat above was closed by
[clone_dosage_depth.md](clone_dosage_depth.md) (job 11784845, 36 samples,
198 clones).

Dosage **is** real — gain-bearing clones are measurably deeper (cell-level
tumor-only rho median +0.104, positive in 24/30 samples, sign test p = 0.0014) —
but it is **far too small** to explain the depth split reported here. Observed
depth spread across tumor clones exceeds the dosage prediction by a median of
**41×**, in **30/30 samples**; `lm(obs ~ pred)` gives slope +10.9, R² = 0.115.
Net dosage range within a sample is a median of 0.020, i.e. a predicted 2% depth
difference, because gains and losses largely cancel.

So the biological alternative I raised against my own design is genuine but
minor. It does not account for the association. It also does not follow that the
association is an artifact — other unmodelled biology (proliferative state, cell
size, total RNA content) remains open. The causal test is still downsampling.
