# Proposed fix: cluster_dictionary over-flagging (excessive cell filtering)

_Status: change #1 (cone/rod threshold 2 → 3) **APPLIED** in numbatHelpers
`R/cluster_dictionary_auto.R`; dictionary regenerated and re-audited. Changes
#2–#4 still proposed. Diagnosis in `doc/cluster_dictionary_audit.md` /
`doc/cluster_dictionary_audit.tsv`._

## Result of the applied cone/rod fix

Raising the cone/rod rule from `>= 2` to `>= 3` of the top-10 markers, then
regenerating the dictionary and re-auditing:

| program | clusters before → after | cells removed before → after |
|---|---|---|
| cone | 23 → **5** | 29,363 → **4,274** |
| rod | 15 → **11** | 6,248 → **2,894** |
| total remove=1 | 64 → **40** | — |
| low-confidence removals | 20 → **3** (only incidental APOE) |

All 5 surviving cone removals now match 3–4 genes (high confidence). The large
cluster-0 deletions on weak `GNB3`+`PCP2`/`RXRG` pairs are gone — e.g. SRX10831286
(was 69.3% removed), SRX10831287 (60.6%), SRX11133593 (64.0%), SRX10264522 (61.0%)
no longer have their main cluster flagged cone. SRX10831280 still loses 82.8% but
that is now genuine MALAT1 (51%) + low_qual/mito (32%) in a 2,991-cell sample, not a
photoreceptor false positive.

## What the audit found

The auto cluster_dictionary flags **64 res.0.2 clusters `remove=1`**, and
`annotate_filter_reason()` drops every cell in those clusters (`cluster_remove`,
highest precedence). This removes an extreme fraction of cells in many samples:

| sample | % cells removed (dictionary only) |
|---|---:|
| SRX10831280 | 82.8% |
| SRX10831286 | 69.3% |
| SRX10264521 | 65.8% |
| SRX10831281 | 65.4% |
| SRX11133593 | 64.0% |
| SRX10264522 | 61.0% |
| SRX10831287 | 60.6% |

Program totals: **cone 23 clusters / 29,363 cells**, low_qual 16 / 18,745,
MALAT1 7 / 11,047, rod 15 / 6,248, APOE 3 / 338.

**Root cause — the `cone` rule fires on 2 weak, non-specific genes.** Of the 23
cone removals, ~15 match the bare 2-gene minimum, and in those the matched genes
are almost always `GNB3` plus `PCP2`/`RXRG` — broadly expressed genes, not
cone-specific. Worse, the flagged cluster is repeatedly **cluster 0, the largest
cluster in the sample**:

| sample | cluster | matched | % of sample | top-10 markers |
|---|---:|---|---:|---|
| SRX10831286 | 0 | GNB3,PCP2 | 69.3% | EYS, OTX2-AS1, GNB3, RCVRN, NEGR1, CKB, FSTL5, KCNB2, PDE6H, SEMA6D |
| SRX10831287 | 0 | PCP2,GNB3 | 59.9% | PCP2, TFF1, LINC02997, EYS, …, GNB3, …, NRL |
| SRX11133593 | 0 | GNB3,RXRG | 56.5% | GNB3, VXN, RCVRN, SNHG29, TFF1, EEF2, RXRG, AIPL1, … |
| SRX10264522 | 0 | GNB3,RXRG | 43.8% | GNB3, VXN, GUK1, CKB, CHCHD2, GSTP1, …, RXRG |
| SRX10264521 | 0 | GNB3,RXRG | 40.7% | VXN, GSTP1, CHCHD2, GUK1, CKB, GNB3, …, RXRG |

These are large, transcriptionally generic clusters — likely tumor/retinal
progenitor populations co-expressing a couple of photoreceptor genes — being
deleted wholesale. By contrast, the 8 high-confidence cone clusters match 3–5
genes from the set, which is the behaviour we want to keep.

`rod` shows the same pattern at smaller scale (e.g. SRX10831287 cl5 matched
`SAG,GNGT1` = 0.68% — a real but tiny rod cluster, fine to keep removed; the
problem is the 2-gene false positives).

`APOE` (3 clusters, 338 cells) is incidental: no rule emits `APOE`, it only
appears via the `top1_gene` catch-all when APOE is rank-1. Those clusters are
actually immune/glia (e.g. SRX11133589 cl4: `APOE, B2M, SPP1, FTL, HLA-DRA,
APOC1, CD74` = macrophage; SRX10264518 cl9: `APOE, VIM, CLU, GPX3, PTN, TF` =
astrocyte). They may belong in the remove set, but the *label/rule* is accidental.

## Proposed changes (numbatHelpers `R/cluster_dictionary_auto.R`)

### 1. Tighten the cone/rod rule (primary fix)

In `classify_cluster_abbreviation()` raise the photoreceptor thresholds from
`>= 2` to `>= 3` of the top-`top_n` markers (lines 65 & 69):

```r
if (length(hit(sets$cone)) >= 3) { ... }   # was >= 2
if (length(hit(sets$rod))  >= 3) { ... }   # was >= 2
```

Per the audit this rescues all ~15 two-gene cone removals (every large cluster-0
deletion above) and the 2-gene rod false positives, while retaining the 8
high-confidence cone and the strong rod clusters (3–8 matched genes).

Optionally also **prune the non-specific genes** from the cone/rod sets
(`.cluster_program_gene_sets`, lines 28–31) — `GNB3`, `PCP2`, `RXRG` drive nearly
every false positive and are not cone-restricted. Dropping them is an alternative
or complement to the threshold change. (Make the threshold change first; it is the
lower-risk, higher-impact lever.)

### 2. Fix the APOE entry

`APOE` in `remove_programs` (line 144) is never produced by a named rule. Either:
- **(a)** drop `"APOE"` from `remove_programs` (its 3 removals are incidental and
  cluster-0-safe at 338 cells), or
- **(b)** add an explicit immune/glia rule (e.g. ≥2 of `APOE, APOC1, CD74,
  HLA-DRA, AIF1, C1QA/B/C`) if these stromal clusters are *meant* to be removed,
  so the label reflects intent rather than a rank-1 accident.

### 3. Tighten the seu-fallback marker call

`.dict_markers_from_seu()` (lines 109–113) ranks by `logFC` after only
`padj < 0.5` — essentially no significance filter, so computed-marker samples get
noisy top genes feeding every rule. Use a real cutoff and a percent-expressed
gate, e.g. `padj < 0.05`, a minimum `logFC`, and `pct_in` above a floor, before
`slice_head(n)`. (Lower priority — most samples use the DB path; this only affects
samples missing from `cluster_markers`.)

### 4. (Consider) the separate cell-level MALAT1 removal

`annotate_filter_reason()` also removes any cell with `abbreviation == "MALAT1"`
independent of `cluster_remove`. With 7 MALAT1 clusters / 11,047 cells this is a
large second lever; revisit `malat1_rank_max` (currently 5) once the cone fix
lands, using the audit to size its impact.

## How to validate a fix once applied (not part of this proposal)

1. Edit the rules above, then re-run `sbatch automate_cluster_dictionary.sbatch`
   to regenerate `results/cluster_dictionary/*.tsv`.
2. Re-run `sbatch audit_cluster_dictionary.sbatch` and diff
   `doc/cluster_dictionary_audit_by_sample.tsv` — the cluster-0 / >40%-removed
   samples should drop sharply; total cone clusters should fall from 23 toward ~8.
3. Only then `tar_invalidate()` the `cluster_dictionary` chain and rebuild
   `filtered_seus` (per CLAUDE.md one-store / seu-cue rules).
