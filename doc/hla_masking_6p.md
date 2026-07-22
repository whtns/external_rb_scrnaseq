# Issue #31: masking HLA genes for 6p inference

**Date:** 2026-07-22
**Verdict: do not run the masking experiment. There is nothing left to mask —
HLA and the entire MHC window are already 100% absent from numbat's 6p
inference, on both the expression and the allele arm. The 6p calls in this
cohort cannot be HLA artifacts.**

Evidence: `results/hla_6p_contribution.csv` (57 non-neutral chr6p segments
across 32 samples), produced by `src/diag_hla_6p_contribution.R`
(job 10488438). No reruns were needed to establish any of this.

---

## 1. The proposed mask would only have reached half the problem

The issue proposes dropping the 27 `HLA-*` genes from the gtf. Tracing how
numbat 1.5.2 actually consumes the gtf, that mask is asymmetric:

- **Expression arm — the mask works.** `run_numbat()` computes
  `genes_annotated = intersect(gtf$gene, rownames(count_mat), rownames(lambdas_ref))`
  and then `count_mat = count_mat[genes_annotated, ]`. Dropping gtf rows
  genuinely removes those genes from the expression model. The issue's finding
  that there is no `blacklist` argument is correct — the gtf *is* the mechanism.
- **Allele arm — the mask does nothing.** `df_allele = annotate_genes(df_allele, gtf)`
  is a `left_join`, so a SNP that no longer hits a gene gets `gene = NA` rather
  than being dropped, and `get_bulk()` then keeps it *deliberately*:
  `filter(lambda_ref != 0 | is.na(gene))`. That `is.na(gene)` clause exists to
  retain intergenic SNPs. Verified empirically: **36.4%** of rows in
  SRX10264524's final bulk already carry `gene = NA` and are retained.

So the mismapping-of-divergent-haplotypes concern — an allele-frequency effect —
would have survived a gtf-only mask untouched. Suppressing it would have needed a
separate positional filter on `allele_df` in `pipeline/scripts/run_numbat.R`.

This turns out to be moot, but it is the reason the experiment as specified would
have produced a misleading "calls survive masking → they're real" result.

## 2. The MHC is already gone — measured, not inferred

Gene survival into the numbat bulk (SRX10031194, representative):

| window | genes in gtf | present in bulk | % |
|---|---|---|---|
| chr6 20.0–28.5 Mb | 134 | 31 | 23.1% |
| **chr6 28.5–33.5 Mb (MHC)** | **248** | **0** | **0.0%** |
| chr6 33.5–40.0 Mb | 94 | 38 | 40.4% |
| genome-wide | 26,807 | 10,674 | 39.8% |

Not one of 248 MHC-window genes survives, and none of the 27 `HLA-*` genes,
in a region whose immediate neighbours retain 23–40%. That is categorical
exclusion, not attrition.

SNP-level, same sample:

| | raw allele counts | final bulk | retained |
|---|---|---|---|
| chr6 20.0–28.5 Mb | 17,001 | 4,986 | 29% |
| **chr6 28.5–33.5 Mb (MHC)** | **16,539** | **9** | **0.05%** |
| chr6 33.5–40.0 Mb | 15,533 | 3,992 | 26% |

The MHC is *denser* than its neighbours in the raw pileup, so this is not a
coverage or phasing-panel gap. And all 9 survivors have `gene = NA` — every SNP
that landed inside an MHC gene was removed. The allele arm is as empty as the
expression arm.

Across all 57 non-neutral chr6p segments in 32 samples, the medians are 0% of
segment expression counts from `HLA-*`, 0% from the whole MHC window, and 0% of
segment SNP depth inside the window.

## 3. Root cause: the reference is cone-only, and cones don't express MHC

`data/sridhar_ref.rds` is a **single-column** reference profile — one cell type,
`Cone`, 15,029 genes. Retina is immune-privileged and photoreceptors are
MHC-low, so the reference assigns essentially zero expression across the MHC.
Two independent numbat filters then remove those genes:

1. `get_bulk()`: `filter(lambda_ref != 0 | is.na(gene))`. 8 of the 18 `HLA-*`
   genes present in the reference have `lambda_ref` exactly 0 and die here.
2. `get_bulk()`: `filter((logFC < 5 & logFC > -5) | Y_obs == 0)`. The remaining
   HLA genes have reference values of 2.3e-07 to 1.7e-04 (`HLA-G` 3.6e-07,
   `HLA-A` 1.4e-04, `HLA-C` 1.7e-04). Any real observed expression against a
   reference that small gives `|logFC| >> 5`. Confirmed active: logFC in the
   bulk is hard-clipped to exactly **[-5.000, 4.990]**.

The allele arm then falls out of the expression arm: once a gene is dropped,
its SNPs still carry the gene *name* but have no `lambda_ref`, so
`filter(lambda_ref != 0 | is.na(gene))` — which spares only `gene = NA` —
removes them too. One filter, both arms.

Note this is *not* the reference simply lacking the genes: 123 of the 248
MHC-window genes and 18 of 27 `HLA-*` genes are present in the reference matrix.
They are present with near-zero values, which is worse than absent, because it
routes them through the logFC filter rather than the membership intersection.

## 4. What this means

- **Issue #31 points 1 and 2 (build a masked gtf, rerun, compare): don't.** The
  comparison would be against an identical result. HLA contributes exactly 0
  genes, 0 counts, and 0 SNP depth to every 6p call already.
- **Points 3 and 4 (widen to the MHC window; 6p-scoped vs cohort-wide): moot**
  for the same reason. Recorded for completeness: the window holds 248 genes of
  which only 27 are `HLA-` prefixed, so a name-prefix mask would have covered
  11% of it — `TRIM27`, `GPX5/6`, `ZBED9`, `TAP1/2`, the olfactory-receptor
  cluster and the `HCG*` lncRNAs would all have been missed. Had masking been
  needed, the positional window was the right unit.
- **For the `SRR13884247` 6p+ review: the 6p+ call is not an HLA artifact.**
  Its 6c segment (8.5–41.1 Mb, LLR 214.2) rests on 101 genes, none of them in
  the MHC.

### The real finding is the flip side

This exclusion is accidental, not principled, and it is not confined to the MHC.
**31.4% of the 15,029 reference genes have `lambda_ref == 0`**, and every one of
them is dropped from inference — plus everything whose tumour expression exceeds
32× the cone value. numbat here only sees genes whose expression sits within a
narrow band of a single cone profile.

The MHC is the most visible casualty because its 5 Mb is contiguous, but the
blind spot is genome-wide, and it systematically removes exactly the genes most
differentially expressed between tumour and normal retina — including
interferon-response and immune genes as a class. Two consequences worth
following up separately:

1. **Segment boundaries drawn across the MHC are interpolated over a ~5 Mb hole
   with zero informative genes and zero informative SNPs.** Several calls bridge
   it: `SRR13884247` 6c (8.5–41.1 Mb), `SRX10031194` 6c (21.2–39.2 Mb),
   `SRR13633762` 6c (21.2–39.2 Mb). The calls are supported by flanking
   evidence; the *breakpoint positions* inside the gap are not.
2. **A multi-cell-type reference would change what numbat can see** cohort-wide,
   not just at 6p. Worth its own issue rather than folding into #31.

## Reproducing

```bash
sbatch diag_hla_6p_contribution.sbatch    # debug partition, ~8 min, read-only
# -> results/hla_6p_contribution.csv, one row per non-neutral chr6p segment,
#    with LLR_x / LLR_y split and HLA / MHC shares of genes, counts, SNP depth
```

The `n_hla` / `n_mhc` columns are 0 throughout. That was checked for a
silent-failure bug and is real: the raw allele counts confirm the MHC SNPs exist
upstream and are removed inside numbat.
