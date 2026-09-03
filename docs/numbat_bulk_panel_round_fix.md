# Round-correct bulk clone panels, with the haplotype track

**Issues:** [#42](https://github.com/whtns/external_rb_scrnaseq/issues/42) (bug),
[#41](https://github.com/whtns/external_rb_scrnaseq/issues/41) (haplotype view)
**Function:** `numbatHelpers::plot_numbat_bulk_clones()`
**Target:** `numbat_bulk_clone_pdfs` → `results/numbat_bulk_clones/<sample>_bulk_clones.pdf`
**Smoke:** `sbatch smoke_numbat_bulk_panel.sbatch`

## The bug

`convert_numbat_pngs()` (`numbat_helpers/R/enrichment_functions_2.R:130`) builds the
bulk-clone panel by globbing `*.png` out of the numbat output directory and
converting `bulk_clones_final.png`. That PNG is numbat's **last** consensus
round — verified byte-identical to `bulk_clones_4.png` (md5 `eaec2970…`,
SRX10264519).

Since `bba4f4a`, `<sample>_numbat.rds` is rebuilt at a **selected** round chosen
by the majority-of-rounds rule, and the heatmap is drawn from that object. So the
bulk panel and the heatmap beside it in every sample summary have been describing
different rounds.

Measured from `results/round_convergence_by_round.csv`:

- **33 / 36** samples: `selected_round != final_round`
- **7 / 36** samples: the canonical RB event set actually differs

| sample | round shown → correct | effect on the panel |
|---|---|---|
| SRX10264524 | 4 → 2 | missing `6p_gain` |
| SRX10264526 | 4 → 1 | missing `13q_loss` |
| SRX11133585 | 4 → 1 | missing `1q_gain`, **spurious `16q_loss`** |
| SRX11133587 | 4 → 2 | **spurious `16q_loss`** |
| SRX11133588 | 4 → 2 | missing `6p_gain` |
| SRX14116946 | 4 → 1 | missing `2p_gain` |
| SRX22868104 | 4 → 3 | missing `16q_loss`, `6p_gain` |

Note the two **spurious** cases: the old panel showed a `16q_loss` that the
selected round does not call.

### Which manifest column to trust

`results/numbat_selected_round.csv` has two columns describing the selected
round's events, and they disagree:

- **`majority_set`** matches the per-round table for **36 / 36** samples. Use this.
- `events_selected` disagrees for **25 / 36**. It appears to carry round 1's
  values. Filed as [#45](https://github.com/whtns/external_rb_scrnaseq/issues/45).

An earlier draft of this work quoted "20 of 39 samples" for the event-set
mismatch. That number came from `events_selected` and was wrong; the correct
figure is 7 of 36.

## The fix

Render from `nb$bulk_clones` — the selected round's own pseudobulk table, present
in **39/39** objects (`results/numbat_rds_qc.csv` never lists `bulk_clones` under
`null_components`) — using numbat's own exported plotter, rather than converting
a PNG off disk.

The consensus round is stamped into the panel title, so this class of mismatch
cannot be invisible again.

## The haplotype page (#41)

The same table already carries the phased-haplotype columns: `pBAF`, `pAD`, `AR`,
`haplo_post`, `haplo_naive`, `haplo_theta_min`, `major_count`, `minor_count`,
`theta_hat`, `theta_mle`, `theta_hat_roll`, `phi_mle`, `p_up`, `loh`.

`numbat::plot_psbulk()` takes `allele_only = TRUE`, so the haplotype view is a
supported argument rather than something we draw ourselves. The function emits
**two pages**: expression + allele, then allele-only (per-clone pHF across the
genome, coloured by CNV state).

This is why #41 and #42 were done together — one function, one data source,
nothing to keep in sync.

Per-cell allelic imbalance (`allele_post_<K>.tsv`: `cell, CHROM, seg, cnv_state,
major, minor, total, MAF, p_loh, …`) would be a separate, heavier panel. Deferred
until the bulk allele track proves useful.

### Gotcha: `gaps_hg38`

`plot_psbulk()` resolves `gaps_hg38` as a **bare symbol** with no argument to
override it. numbat ships it as `LazyData`, which lives in the package's
lazy-load database and is reachable only once numbat is **attached** — it is not
in the namespace, so importing is not enough. Without attaching, every panel dies
with `object 'gaps_hg38' not found`. The function attaches numbat once,
idempotently.

## The duplicated targets

`large_numbat_pdfs`, `filtered_numbat_pdfs` and `low_hypoxia_numbat_pdfs` had
**identical commands** — all three were `convert_numbat_pngs(numbat_rds_files)` —
so the work ran three times per sample to produce the same files.

There is no filtered or low-hypoxia numbat run for any in-scope sample:
`output/numbat_sridhar_filtered/` holds 27 directories, **all SRR, zero SRX**,
and there is no low-hypoxia numbat directory at all. Building them would mean
running numbat 39×2 more times, which is out of scope (and forbidden for
`SRX11133592/93/94`).

**Decision: keep duplicating the unfiltered run, but stop computing it three
times and label it honestly.** One target now feeds all three columns, and rows 4
and 5 are labelled *"(numbat run; same in all columns)"*.

The three-column layout is unchanged, because the other slots genuinely differ
per cell set:

| slot | source | differs per column? |
|---|---|---|
| count strip, karyogram, clone/segment trees | per-cell-set | yes |
| heatmaps | `numbat_heatmap_plots_{unfiltered,subset,low_hypoxia}`, built from `unfiltered_seus` / `filtered_seus` / `seus_low_hypoxia` | yes |
| expression, bulk clones | one numbat run per sample | **no — deliberately duplicated** |

## Verification

Smoke test on SRX10264524, SRX22868104 (two of the seven event mismatches) and
SRX10264519 (control):

```
   sample_id pdf_ok n_pages round_held round_selected round_final
 SRX10264524   TRUE       2          2              2           4
 SRX22868104   TRUE       2          3              3           4
 SRX10264519   TRUE       2          2              2           4
```

All three render two pages, and the round each object holds equals
`selected_round`. Object events match `majority_set` in every case.

Rendering costs ~20–47 s per sample, which is why the target runs on the heavy
crew controller.
