# Issue #31: masking HLA genes for 6p inference

**Date:** 2026-07-22
**Verdict: numbat already masks the MHC, by design. `numbat:::filter_genes()`
hardcodes an exclusion of chr6:28,510,120-33,480,577, which removes exactly 248
genes including all 27 `HLA-*`. There is no masking work to do, and the 6p calls
in this cohort cannot be HLA artifacts.**

Evidence: `results/hla_6p_contribution.csv` (57 non-neutral chr6p segments across
32 samples), from `src/diag_hla_6p_contribution.R` (job 10488438). No reruns
were needed.

> **Revision note.** An earlier version of this document attributed the missing
> MHC to the reference profile being a single-cell-type (`Cone`) matrix. That was
> wrong, and is corrected in section 3. It came from reading
> `data/sridhar_ref.rds` (15,029 x 1, Cone-only) — which **no run actually
> used**. The runs used
> `/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds`, a different file
> with the same name: 22,901 x 24, twenty-four retinal cell types. The
> measurements in sections 1-2 were direct observations of numbat output and are
> unaffected.

---

## 1. numbat's default gene blacklist is the MHC

The issue concluded numbat 1.5.2 has no default blacklist, having searched for a
`blacklist` / `exclude_genes` argument and checked the exported data objects.
Both checks are correct as far as they go — but the exclusion is neither an
argument nor exported data. It is a hardcoded coordinate filter inside the
unexported `numbat:::filter_genes()`, called from `get_exp_bulk()`:

```r
genes_exclude = gtf %>%
  filter(CHROM == 6 & gene_start < 33480577 & gene_end > 28510120) %>%
  pull(gene)
genes_keep = genes_keep[!genes_keep %in% genes_exclude]
```

Applied to `gtf_hg38` this drops **248 genes, including 27 of 27 `HLA-*`** — not
one HLA gene escapes it. numbat's authors evidently reached the same conclusion
the issue does, and acted on it at the region level rather than the gene-name
level, which is what the issue's point 3 argues for.

## 2. Confirmed in the output: the MHC contributes nothing to any 6p call

Gene survival into the numbat bulk (SRX10031194, representative):

| window | genes in gtf | in bulk | % |
|---|---|---|---|
| chr6 20.0-28.5 Mb | 134 | 31 | 23.1% |
| **chr6 28.5-33.5 Mb (MHC)** | **248** | **0** | **0.0%** |
| chr6 33.5-40.0 Mb | 94 | 38 | 40.4% |
| genome-wide | 26,807 | 10,674 | 39.8% |

The gap is surgical: the bulk runs to `ZSCAN12` at 28.38 Mb, then jumps to
`BAK1` at 33.57 Mb, with flanking genes surviving at ordinary logFC (-2.1 to
+4.1). That is the blacklist boundary, not attrition.

The allele arm goes with it. SNPs inside a blacklisted gene keep the gene *name*
from `annotate_genes()` but have no `lambda_ref`, so `get_bulk()`'s
`filter(lambda_ref != 0 | is.na(gene))` — which spares only `gene = NA` —
removes them:

| | raw allele counts | final bulk | retained |
|---|---|---|---|
| chr6 20.0-28.5 Mb | 17,001 | 4,986 | 29% |
| **chr6 28.5-33.5 Mb (MHC)** | **16,539** | **9** | **0.05%** |
| chr6 33.5-40.0 Mb | 15,533 | 3,992 | 26% |

All 9 survivors have `gene = NA`, i.e. they fall in intergenic gaps. The MHC is
*denser* than its neighbours in the raw pileup, so this is not a coverage or
phasing-panel gap.

Across all 57 non-neutral chr6p segments in 32 samples: median 0% of segment
expression counts from `HLA-*`, 0% from the whole MHC window, 0% of segment SNP
depth inside it.

**For the `SRR13884247` 6p+ review (#30): the call is not an HLA artifact.** Its
6c segment (8.5-41.1 Mb, LLR 214.2) rests on 101 genes, none in the MHC.

## 3. The reference is fine — correcting the earlier claim

The reference actually used,
`/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds`, is **22,901 x 24**
with 24 retinal cell types (including Macrophage, Microglia, Monocyte, T/NK-Cell,
which do express MHC). 7.8% of its genes are zero across all columns.

`choose_ref_cor()` then picks one column per cell, and assigns the overwhelming
majority to `Cone` — 1397/1458 for SRX10031194, 9031/10000 for SRR13884247,
10611/11643 for SRX10264524 (90-96%). That is expected for a cone-precursor
derived tumour, not a defect. HLA class I is low but **nonzero** in that column
(`HLA-A` 48.3e-6, `HLA-C` 56.2e-6, `HLA-B` 7.5e-6), and only 28.8% of
MHC-window genes are zero there — so the reference does **not** explain the
0/248 wipeout. The hardcoded blacklist does, completely and exactly.

**Trap worth knowing:** `data/sridhar_ref.rds` and `data/sridhar_normal_ref.rds`
are byte-identical (md5 `7f34cbed…`), both 15,029 x 1 Cone-only, both dated
2022-10-14, and **neither is used by any run**. Every `done.txt` in
`output/numbat_sridhar/` records either the `/project2/.../Homo_sapiens/numbat/`
path (37 samples) or `/dataVolume/storage/...` from the original machine (32
samples). Reading the `data/` copy to find out what the runs used gives the wrong
answer.

## 4. What this means

- **Issue #31 points 1-2 (build a masked gtf, rerun, compare): don't.** The mask
  is already applied upstream of everything, by numbat itself.
- **Point 3 (whole MHC window vs `HLA-` prefix): correct, and already how numbat
  does it.** Recorded for completeness: the window holds 248 genes of which only
  27 are `HLA-` prefixed, so a name-prefix mask would have covered 11% —
  `TRIM27`, `TAP1/2`, `GPX5/6`, `ZBED9`, the olfactory-receptor cluster and the
  `HCG*` lncRNAs all missed. numbat excludes the region, catching all of them.
- **Point 4 (6p-scoped vs cohort-wide): cohort-wide, already.**

### One residual, real consequence

Segment boundaries drawn across the blacklisted window are interpolated over
~5 Mb containing zero informative genes and zero informative SNPs. Three calls
bridge it: `SRR13884247` 6c (8.5-41.1 Mb), `SRX10031194` 6c (21.2-39.2 Mb),
`SRR13633762` 6c (21.2-39.2 Mb). The calls themselves are supported by flanking
evidence and stand; the **breakpoint positions inside the gap are not
determined** and should not be quoted to Mb precision. This is a property of
numbat's deliberate design choice, not a bug to fix here.

### Method note, for the record

Had masking *not* already been done, the gtf-only mask in point 1 would have been
asymmetric: it reaches the expression arm (`count_mat` is subset to
`intersect(gtf$gene, ...)`), but not the allele arm, because `annotate_genes()`
is a `left_join` and `get_bulk()` deliberately retains `gene = NA` SNPs — 36.4%
of SRX10264524's bulk rows already are. Removing a gene from the gtf converts its
SNPs into retained intergenic SNPs rather than dropping them. The experiment
would have returned "calls survive masking -> they're real" for the wrong reason.
Suppressing the allele arm needs a positional filter on `allele_df` in
`pipeline/scripts/run_numbat.R`.

## Reproducing

```bash
sbatch diag_hla_6p_contribution.sbatch    # debug partition, ~8 min, read-only
# -> results/hla_6p_contribution.csv
```

The `n_hla` / `n_mhc` columns are 0 throughout. That was checked for a
silent-failure bug and is real: the raw allele counts confirm the MHC SNPs exist
upstream and are removed inside numbat.
