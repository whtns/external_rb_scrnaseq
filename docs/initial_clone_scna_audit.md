# Which tumors need the aggregate diploid object as their comparator?

**30 of 39 samples have at least one canonical RB SCNA already present in their
earliest tumor clone.** In those samples the SCNA is clonal — every tumor cell
carries it — so there is no within-sample SCNA-negative tumor clone to compare
against, and its effect cannot be measured internally. The aggregate diploid
object is the only available comparator. This is the majority of the cohort, not
an edge case.

| RB SCNAs already in the founder tumor clone | samples |
|---|---|
| 3 | 6 |
| 2 | 13 |
| 1 | 11 |
| 0 — founder is RB-SCNA-free | 9 |

Per-sample table: `results/diploid_audit/founder_clone_rb_scna.csv`, with
`rb_in_founder` (clonal — needs the aggregate) and `rb_acquired_later`
(subclonal — internally testable) split out for each sample. Those two columns
are the practical output: SRX11133592 and SRX11133594, for instance, carry `16q-`
in the founder but acquire `1q+` later, so `1q+` can be tested within-sample there
while `16q-` cannot.

### A note on how this was measured

An earlier version of this document reported "0 of 39 samples have an initial
clone carrying an SCNA". That measured the wrong thing. numbat's
`clone_opt == 1` is *always* the empty-`GT_opt` normal clone — it is the model's
reference by construction, not a biological finding — so that question can only
ever answer zero. The number above instead takes the **earliest tumor clone**,
the clone with the fewest `GT_opt` tokens, which is what "initial clone with an
SCNA" means in practice.

### The in-sample comparator is usually too small anyway

Even where numbat did find a normal clone, it is rarely big enough to serve as a
differential-expression comparator. Among the 22 panel samples:

- **9 have fewer than 50** diploid cells
- **14 have fewer than 200**

with SRX14116947 at 2 cells and SRX22868104 at 6, against 1,073 and 12,745
low-hypoxia cells respectively. So the aggregate diploid is needed nearly
cohort-wide, not only for the 30 samples with a clonal RB SCNA.

### Why this raises the stakes on the contamination finding below

If the aggregate object is the comparator for 30 samples, its purity is no longer
a housekeeping concern — tumor cells pooled into it bias every one of those
comparisons toward the null. The rest of this document is about exactly that, and
about the fix that makes the object usable for this purpose.

### One caveat on the comparison itself

Comparing a sample's tumor cells against diploid cones pooled from *other*
samples confounds the SCNA effect with donor and batch. `assemble_diploid_seu()`
integrates across samples (`seuratTools::integration_workflow`), which mitigates
this but does not remove it. A clonal-SCNA effect measured this way is not the
same quantity as a within-sample clone contrast, and the two should not be
pooled uncritically.

## But the audit surfaced a different way tumor cells reach the diploid pool

`assemble_diploid_seu()` (`plot_functions_3.R:653`) selects its diploid cells as

```r
diploid_mask <- is.na(seu$scna) | seu$scna == ""
```

and `seu$scna` is derived per cell from that cell's `GT_opt`
(`plot_functions_3.R:150-171`). So a cell reaches the pool by three routes, and
only the first means "diploid":

| route | how | intended? |
|---|---|---|
| **A** | the cell's clone has an empty `GT_opt` — numbat's normal clone | yes |
| **B** | the cell is absent from `clone_post`, so `scna` is `NA` | no, silent |
| **C** | the clone has a **non-empty** `GT_opt` whose tokens do not resolve against the clone-simplification key | no, silent, and aneuploid |

Route B is clean: **0 cells** cohort-wide. Every cell in every
`*_hypoxia_low_seu.rds` is present in its sample's `clone_post`.

Route C is not. `simplify_gt_col()` (`plot_functions_15.R:236-250`) joins the
clone's `GT_opt` tokens against the key, drops tokens that do not match, and
then:

```r
if (length(gt_vals) == 0) gt_vals <- ""
```

An unresolved clone collapses to `""` — **the same value a genuinely diploid
clone gets**. There is no distinct "unknown" state, so a fully aneuploid clone
becomes indistinguishable from a normal one downstream.

### Why the founder clone is the one that loses

In all affected samples the mislabelled clone is the **founder tumor clone** —
the clone with the fewest events. That is systematic, not chance.
`compute_clone_simplifications()` (`metadata_functions_1.R:507-509`) keeps **one
representative segment per SCNA label**, the alphabetically first:

```r
result <- non_neu[order(non_neu$scna_label, non_neu$seg), ]
result <- result[!duplicated(result$scna_label), ]
```

A label like `16q-` may be carried by several segments (`16b`, `16d`, `16e`,
`16f`); only one enters the key. Descendant clones accumulate enough segments to
hit the representative and resolve; the founder, holding only its own initial
segments, often does not.

SRX22868103 at its selected round shows it exactly (`n_tokens_resolved` from the
per-clone table):

| clone | `GT_opt` | events | cells | tokens resolved |
|---|---|---|---|---|
| 1 | *(empty)* | — | 607 | 0 — genuinely diploid |
| **2** | `16e,1q,1n,1s` | 16q-;1q+ | **596** | **0 → labelled `""`** |
| 3 | `16e,1q,1n,1s,11b,6a` | 11p-;16q-;1q+;6p+ | 896 | 2 |
| 4 | `16e,1q,1n,1s,16d,16f` | 16q-;1q+ | 3822 | 1 |
| … | | | | |

Clone 2 carries `16q-` and `1q+` — the two canonical RB events — and is labelled
diploid.

A second, compounding detail: manual entries in
`config/large_clone_simplifications.yaml` are merged with
`modifyList(computed, manual)` (`R/pipeline_targets_inputs.R:236`), which
**replaces** the computed representative for a label rather than adding to it.
For SRX22868103 the YAML pins `16q-: 16b` and `1q+: 1h`, neither of which clone 2
carries. Curating the YAML can therefore *narrow* coverage. Auditing with the
computed key alone finds 3 affected samples; with the merged key the pipeline
actually uses, it finds 6.

### Scale

Route-C cells are counted only within each sample's `hypoxia_low` cell set — the
input `diploid_seu` actually draws from — before the cone-type restriction.

| sample | in panel | pool A (true diploid) | pool C (mislabelled) | % contaminated |
|---|---|---|---|---|
| SRX10264525 | TRUE | 39 | 9094 | 99.6% |
| SRX22868103 | TRUE | 210 | 4107 | 95.1% |
| SRX11133589 | TRUE | 1627 | 584 | 26.4% |
| all other panel samples | TRUE | 4133 | 0 | 0% |

**Cohort pool: 6,009 genuinely diploid cells vs 13,785 aneuploid ones — 70% of
the pool would be tumor cells,** concentrated in two samples that would each
contribute more mislabelled tumor cells than the entire rest of the cohort
contributes real ones.

Three further samples are affected but do not reach the pool: SRX11133585 and
SRX14116944 are `include = FALSE`, and SRX11133590 is not in the panel table.

## Also found: two samples silently dropped from the panel

`diploid_panel_sample_ids` filters with `%in%`, so a sample present in
`seus_low_hypoxia` but absent from `data/diploid_panel_samples.csv` is excluded
with no recorded reason — the exact failure mode #15 asks to eliminate. Two
samples are in that state:

- **SRX11133587**
- **SRX22868105**

Both should get an explicit row with a reason, whatever the decision.

## The blocker: `scna` is empty for every cell in every sample

The route-C counts above are the floor. The ceiling is worse, and it is what
actually holds.

`diag_scna_label_reality.sbatch` (job 11785042) loaded three saved
`*_hypoxia_low_seu.rds` and tabulated `scna` directly:

| sample | cells | distinct `scna` values | cells with `GT_opt` carrying an SCNA but `scna == ""` |
|---|---|---|---|
| SRX22868103 | 12,472 | `""` only | 12,262 |
| SRX10264525 | 9,133 | `""` only | 9,094 |
| SRX10264519 | 3,630 | `""` only | 3,527 |

**Every cell in every object is labelled `""`.** SRX10264519, which the
TSV-level audit scored as completely clean, is affected too. `assemble_diploid_seu()`'s
mask matches 12,472 of 12,472 cells: it is a no-op, and the "consensus diploid"
object built from it would be the entire low-hypoxia cohort.

### Root cause

`large_clone_simplifications_per_sample` (`R/pipeline_targets_inputs.R:231-241`)
returns the merged key directly:

```r
if (is.null(manual) || length(manual) == 0) computed else modifyList(computed, manual)
```

That value is a **`label -> seg` list** (`compute_clone_simplifications()` sets
`names(out) <- result$scna_label`). Both consumers then index it **by sample id**:

```r
if (is.null(large_clone_simplifications[[sample_id]])) {          # plot_functions_3.R:142, :278
  large_clone_simplifications <- tibble::tibble(scna = character(), seg = character())
```

`"SRX10264519"` is not a label, so the lookup returns `NULL`, the key falls back
to the empty tibble, every token drops, and every clone collapses to `""`.

The sibling target gets this right — `large_filter_expressions_per_sample`
(`:246`) slices with single brackets, `large_filter_expressions[sample_id]`,
which keeps the sample name so `[[sample_id]]` resolves. The function's own
docstring states the intended contract: *"Named list of per-sample SCNA
simplification vectors"*. The target is what is wrong, not the functions.

The merge target dates from `14d58c0` (2026-06-05); the objects are from Aug 26,
so they were built with this defect. Both `unfiltered_seus` and `filtered_seus`
are affected identically.

## What was done about it

**The diploid object no longer depends on any of this.** `assemble_diploid_seu()`
now selects on `GT_opt`, numbat's own genotype string, which is stored on the
same objects and is empty only for numbat's normal clone:

```r
diploid_mask <- !is.na(seu$GT_opt) & seu$GT_opt == ""
```

`NA` (cell absent from `clone_post`) is treated as unknown and excluded rather
than assumed normal — currently 0 cells cohort-wide, and the guard keeps it so.
A missing `GT_opt` column is a hard error, so the function can never silently
fall back to a label that means two different things.

This makes the object correct **without rebuilding `filtered_seus`**, because
`GT_opt` is already present on the saved objects.

## Two things left for you to decide

**1. The key-shape bug is not fixed.** The one-line correction is to return a
list keyed by sample id, matching the sibling target:

```r
setNames(list(merged), sample_id)
```

I have not applied it. Changing that target's value invalidates
`unfiltered_seus` and `filtered_seus` cohort-wide, which is a full-pipeline
rebuild — that is your call to schedule, not a side effect I should arm silently.
Until it lands, `seu$scna` is `""` everywhere, so **every figure keyed on `scna`
is currently showing one undifferentiated clone**. Worth checking how far that
reaches.

**2. Even with that fixed, route C survives for 6 samples.**
`compute_clone_simplifications()` keeps one representative `seg` per label, and
`GT_opt` tokens are `seg_cons` values, not `seg`. A founder clone whose segments
are not the chosen representatives still resolves zero tokens and still collapses
to `""`. Fixing it properly means keying the map on `seg_cons` and retaining
every segment per label (deduplicating labels when several match), and making
`simplify_gt_col()` emit something other than `""` when nothing resolves, so
"unlabelled" is never confusable with "diploid". That changes label text in
existing figures, which is why it is a decision rather than a patch.

## What this audit does not establish

It does not judge whether the diploid cells are cones —
`assemble_diploid_seu()` restricts to cone cells after the diploid mask, so the
final object is smaller than the pool counts here. Nor does it revisit which
samples belong in the panel; it only reports that two are missing from the table.
