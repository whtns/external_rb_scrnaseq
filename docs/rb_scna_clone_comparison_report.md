# Which tumors permit a with/without RB SCNA clone comparison?

Report page: https://claude.ai/code/artifact/d63c5250-3273-4827-9a25-075ca3af4808

**16 samples have a clone pair differing by an RB SCNA while both clones stay
aneuploid** (25 distinct sample x SCNA contrasts). Eight of the nine samples of
interest are confirmed; `SRX10831287` is not.

## Tiers

- **Tier A** — both clones aneuploid, i.e. the SCNA is subclonal. The tumor
  background is held fixed, so the contrast isolates the event. This is the
  comparison of interest.
- **Tier B** — the only SCNA-negative clone is numbat's diploid clone, i.e. the
  SCNA is clonal. The contrast is tumor-vs-normal and confounds the event with
  everything else separating tumor from normal cells. These samples need the
  aggregate diploid object instead (see `initial_clone_scna_audit.md`).

Eligibility is decided by nesting: one clone's `GT_opt` token set is a subset of
another's, and the added tokens map to an RB SCNA. Both clones must clear 20
cells. Computed at the **selected** round from `clone_post` / `segs_consensus`.

## Correction to the first version

The first run tested nesting on `seg_cons` tokens alone, then mapped the gained
tokens to arm labels. That counts a clone gaining a **second segment on an arm it
already carries** as acquiring the SCNA. `SRX22868103` clones 2 and 4 are both
`16q-;1q+`; `SRX11133593` clones 2 and 3 are both `16q-;1q+`; `SRX10264526`
clones 2 and 3 are both `16q-;1q+;2p+`. Nested at token level, but not
with/without contrasts.

Eligibility now requires the **arm-level label to be absent from the preceding
clone and present in the descendant**. The set drops from 16 samples / 25
contrasts to **11 samples / 15 contrasts**, and three of the nine listed samples
are removed.

## Answers to the question asked

| | |
|---|---|
| Confirmed tier A, of the nine | SRX10264523 (2p+), SRX11133592 (1q+), SRX11133594 (1q+), SRX22868103 (6p+), SRX14116944 (6p+), SRX14116947 (1q+) |
| Not eligible | **SRX10264526, SRX10831287, SRX11133593** |
| Additional tier A samples | SRX10031193, SRX10831282, SRX11133588, SRX14116946, SRX22868105 |

**Note the SCNA, not just the sample.** Of the seven listed under "1q and 16",
none has a testable 16q- contrast: 16q- is clonal in all of them. The testable
events are 1q+ (SRX11133592, SRX11133594), 2p+ (SRX10264523) and 6p+
(SRX22868103). Of the two listed under "2p or 6p", SRX14116944 tests 6p+ but
SRX14116947 tests only 1q+ — both its 2p+ and 6p+ are clonal.

## Tier A contrasts (15)

| SCNA | samples | clean (no co-event) |
|---|---|---|
| 1q+ | 6 | 3 |
| 6p+ | 5 | 1 |
| 16q- | 2 | 0 |
| 2p+ | 2 | 0 |

16q- is the scarcest: clonal in 14 samples, with only `SRX22868105` (3,422 vs
1,855 cells) well powered internally. Full table:
`results/diploid_audit/rb_scna_comparison_eligibility.csv`.

## Tier B — fit for the aggregate diploid object

26 samples carry at least one RB SCNA clonally (1q+ in 16, 16q- in 14, 2p+ in 7,
6p+ in 7), 44 sample x SCNA contrasts with no internal comparator. Largest:
SRX10264525 (6p+, 12,583 cells), SRX10831287 (1q+, 8,452), SRX14116945
(16q-;1q+, 6,190), SRX10264521 (1q+, 4,283), SRX22868103 (16q-;1q+, 3,822).

**The aggregate diploid object is not yet fit for this.** It holds 3,117 cone
cells from 16 samples, but SRX11133589 alone is 51% of it, and two large
contributors (SRX10264526, 835 diploid cells; SRX11133594, 653) dropped to zero
because none of their diploid cells typed as cones.

## Config changes applied

Five comparison keys added to `config/large_clone_comparisons.yaml` (backup:
`.bak_20260907`); no existing key changed:

| Sample | Key | Segments |
|---|---|---|
| SRX10031193 | `3_v_2_1q+` | 1i, 1n, 1k |
| SRX10031193 | `3_v_2_2p+` | 2c, 2f |
| SRX10031193 | `3_v_2_6p+` | 6c, 6b |
| SRX11133588 | `3_v_2_16q-` | 16c |
| SRX11133588 | `3_v_2_1q+` | 1g |

`scna_collage_samples` extended: 1q 6->10, 6p 8->9, 16q 3->5, 2p unchanged at 9.
Existing tier-B entries kept — the tumor-vs-diploid collage is still wanted, it
is simply a different contrast and should be read as such.

## Artifact status

- **Sample summaries**: present for all tier-A samples except SRX10031193 and
  SRX11133591.
- **Two-clone collages**: four merged per-SCNA PDFs, built Aug 28 against the
  OLD sample sets; they do not yet include the added samples.
- **SCNA-specific diffex**: only 2 of 15 tier-A contrasts, both 1q+
  (SRX11133592, SRX11133594). None for 16q-, 2p+ or 6p+.

## Caveats

Nesting is judged at the selected round. Cell counts are from `clone_post`,
pre-filtering. The 20-cell floor admits pairs too thin for reliable diffex —
SRX11133588 (29 / 66 cells) is flagged, not excluded.

Tables: `results/diploid_audit/rb_scna_comparison_eligibility.csv`,
`clone_pairs_rb_scna.csv`. Scripts: `src/audit_clone_pairs_rb_scna.R`,
`src/gen_tierA_clone_comparisons.R`, `src/build_rb_scna_comparison_report.R`.
