# Maximizing numbat recovery of RB SCNAs — parameter recommendations

**Date:** 2026-08-19 (updated 2026-08-20 with the pilot result)

> **Bottom line after the pilot.** The only recommendation here that survived
> testing is §2, the cross-round union — it recovers **+37 canonical events
> (172 → 209, +22%)** across the cohort at zero compute cost. `multi_allelic =
> TRUE`, originally billed as the highest-yield rerun, was tested and **loses**
> events (§3a); leave it `FALSE`. The `t` question is closed for lack of
> same-platform data (§4).
**Scope:** the 71-sample production cohort in `output/numbat_sridhar/`, numbat 1.5.2
**Basis:** prior parameter tests in `doc/diploid_baseline_pin_6p.md`,
`doc/hla_masking_6p.md`, `src/diag_consensus_convergence.R`,
`src/sweep_6p_gain_candidates.R`, plus two new diagnostics written here:
`src/diag_rb_scna_recovery.R` and `src/diag_numbat_param_regimes.R`.

Canonical RB arm events throughout are 1q gain, 2p (MYCN) gain, 6p gain,
13q loss, 16q loss. hg38 boundaries and state matching are documented at the top
of `src/diag_rb_scna_recovery.R`.

---

## 1. Two constraints that set the scope

**Only 39 of 71 samples can be rerun at all.** The 32 SRR samples were run on the
previous machine and their inputs — cellranger matrix, seurat object, allele
counts — were never migrated here. All three are absent for every one of them.
Rerunning those means regenerating from FASTQ. Every rerun-based recommendation
below therefore covers the 39 SRX samples only; the cross-round union in §2 is
the sole lever that reaches the whole cohort.

**The cohort was not run under one parameter set.** numbat echoes its real
parameters into `output/numbat_sridhar/<sample>/log.txt`, and those disagree with
`pipeline/config.yaml`, which has drifted. Actual regimes
(`results/numbat_param_regimes.csv`):

| t | alpha | max_iter | min_LLR | max_entropy | n | of which SRX |
|---|---|---|---|---|---|---|
| 1e-2 | 1e-3 | 4 | 5 | 0.7 | 31 | 31 |
| 1e-2 | 1e-4 | 4 | 5 | 0.7 | 6 | 6 |
| 1e-2 | 1e-4 | 2 | 5 | 0.5/0.7/0.8 | 16 | 0 |
| 1e-2 | 1e-4 | 1 | 5 | 0.5 | 1 | 0 |
| 1e-5 | 1e-4 | 2 | 2 | 0.5/0.8 | 17 | 2 |

`multi_allelic = FALSE` in all 71. Config currently says `alpha: 1e-4`; 31 of the
39 rerunnable samples actually ran at `1e-3`. **Any rerun must replay each
sample's logged parameters, not config.yaml**, or it confounds the change under
test with an unintended alpha change.

Usefully, the 39 rerunnable samples are nearly homogeneous: 37 at
`t = 1e-2, max_iter = 4`, 2 at `t = 1e-5`.

## 2. The free recovery: score segments across consensus rounds

numbat's iteration does not converge here — only 4/71 runs reach `d_k = 0`, and
freezing the diploid baseline did not change that
(`doc/diploid_baseline_pin_6p.md`). So the final `segs_consensus_<K>.tsv` is a
stopping point, not an answer, and a call present in rounds 1–3 and absent at
round 4 is not evidence against the call.

Counting canonical events across all 71 samples (`src/diag_rb_scna_recovery.R`):

| | events |
|---|---|
| final consensus round (what the pipeline consumes) | **172** |
| union across all rounds of the same runs | **209** |

**+37 events (22%), already on disk, no compute.** Per event (all 71 dirs, before
the provenance exclusion below):

| event | final | union | recovered, stable | recovered, marginal |
|---|---|---|---|---|
| 1q gain  | 48 | 53 | 1 | 4 |
| 2p gain  | 26 | 36 | 1 | 9 |
| 6p gain  | 25 | 35 | 2 | 8 |
| 13q loss | 34 | 44 | 4 | 6 |
| 16q loss | 43 | 50 | 2 | 5 |

Read the tiering honestly: of the 37, only **9 are stable** (called in ≥50% of
that sample's rounds); **28 are single-round** calls. This is a recovery of 9
solid events plus 28 candidates that need a second line of evidence — not 37 free
confident calls. The worst single case is `SRX10264524`: one event at round 4,
all five across rounds 1–4.

**Provenance guard — why 172/209 and not 176/218.** The union assumes
`segs_consensus_1..K` in a sample dir are successive rounds of *one* run. In four
dirs they are not: a later rerun crashed partway and overwrote only the early
rounds, leaving the older run's later rounds in place. `SRX10031192` round 1 is
dated 2026-06-03 and round 2 is 2026-04-20; `SRR13884242` round 1 is 2024-03-18
against rounds 2–3 on 2024-02-19. In those dirs a "round-to-round" difference is
partly a run-to-run difference under possibly different parameters, so the
stability score is not interpretable. `src/diag_rb_scna_recovery.R` now flags any
dir whose round mtimes span more than a day (`mixed_provenance` column) and
reports clean dirs separately: **SRR13884242, SRR27187901, SRX10031192,
SRX11133590**, 9 events between them. Quote the clean figure.

**A related but separate defect: stale `log.txt`.** Three more dirs
(`SRX10031191`, `SRR17960483`, `SRR27187899`) hold a `log.txt` from a *later* run
than their segments — for `SRX10031191` the log is from 2026-06-04 and stops at
"Iteration 1" while the segments are from 2026-03-24. The parameters in those
logs describe a run that produced nothing, so `src/diag_numbat_param_regimes.R`
mis-attributes their regime and `rerun_numbat_arm.sbatch` would replay a crashed
run's settings. Check `ITERS` in the log against the number of
`segs_consensus_*.tsv` before trusting a dir's parameters.

**Design constraint on how to consume this.** The union cannot be retrofitted
into the `<sample>_numbat.rds` that the targets pipeline reads.
`pipeline/scripts/process_numbat_rds.R` packages `joint_post`, `exp_post`, and
`allele_post` at a single `target_iter`, and those are cell-level posteriors keyed
on that round's `seg` labels — which numbat re-letters every round. A union
segment table would no longer align with them. Keep it as a scoring and reporting
layer beside the rds, not a replacement inside it.

```bash
Rscript src/diag_rb_scna_recovery.R
# -> results/rb_scna_segments_by_round.csv    every non-neutral segment,
#                                             n_rounds_called / n_rounds, max LLR
#    results/rb_scna_canonical_by_sample.csv  sample x event, final vs union, tier
```

## 3. `multi_allelic = TRUE` — tested and rejected

> **Mechanism claim retracted 2026-08-20.** This section originally argued that
> `multi_allelic = FALSE` disabled `test_multi_allelic()` so no segment could
> resolve to `bamp`/`bdel`. That is false: **24 of 71 production dirs already
> contain `bamp`/`bdel` at their final round** under `multi_allelic = FALSE`, and
> on `SRX10031192` numbat's own log under `TRUE` reports `0 multi-allelic CNVs
> found`. Whatever `multi_allelic` does here, it is not unlocking otherwise
> unreachable states. The pilot result below supersedes the argument; see §3a.

Production hardcoded `FALSE` at the `run_numbat()` call site. Validated on two
samples at matched round 4:

| sample | production | `multi_allelic = TRUE` |
|---|---|---|
| `SRX10264524` | 1 canonical event (13q loss) | **4** (+6p, +1q, +2p) |
| `SRX14116944` | 3 (6p, 16q, 1q) | **4** (+2p) |

chr6p on `SRX10264524` goes from absent to LLR ~300 at 21.2–41.1 Mb; chr1q from
LLR 65.0 → 433.0 on identical coordinates. Not uniform — `SRX14116944` chr6p
0.3–49.5 Mb drops 524.9 → 465.9 on identical coordinates — so the case rests on
recovering calls lost outright, not on every LLR rising.

`pipeline/scripts/run_numbat.R:41-42` already defaults to TRUE and the Snakefile
does not override it, so this needs no code change.

**Both rows of that table are invalid.** They compare `multi_allelic = TRUE`
against production's *final round*, not against production's *cross-round union*
(§2) — which costs nothing and is already available. Scored against the union:

| sample | production union | `multi_allelic = TRUE` | real delta |
|---|---|---|---|
| `SRX10264524` | **5** (13q, 16q, 1q, 2p, 6p) | 4 | **-1** |
| `SRX14116944` | **4** (16q, 1q, 6p, 2p) | 4 | **0** |

`SRX10264524` calls all five canonical events somewhere in its four production
rounds. The apparent "+3" was the union recovering them, not `multi_allelic`.

## 3a. Pilot result — `multi_allelic = TRUE` loses events

Six samples (alphabetical, `SRX11133592/93/94` excluded), replaying each sample's
own logged parameters and changing only `multi_allelic = TRUE` and `max_iter = 4`.
Job `11228996`, `output/numbat_ma/`. Task 0 (`SRX10031191`) failed for an unrelated
reason (§3b); five samples scored.

| | final | **union** |
|---|---|---|
| production | 8 | **10** |
| `multi_allelic = TRUE` | 9 | **9** |

Zero canonical events gained. One lost.

- **`SRX10264518` 6p gain — lost.** Production calls a focal `amp` on 6p in
  rounds 1 (27–65 Mb, LLR 95.9) and 3 (24–43 Mb, LLR 62.4). Under `TRUE` chr6
  collapses from 6–8 segments per round to 3–4, and everything from 65 Mb to
  170 Mb — across the centromere — merges into one `loh`/`del` block whose LLR
  runs away across rounds: 287.8 → 428.7 → 1513.3 → 1584.8. The focal 6p gain is
  absorbed into that block and never called.
- **`SRX10031192` 6p gain — cosmetic.** It moves from union-only to final-round,
  but that sample halts after one iteration ("No CNV remains after filtering by
  entropy in single cells"), so final and union are the same round by definition.

The mechanism is over-merging, not state unlocking. Across the five samples'
final rounds, `bamp`/`bdel` barely move (7/1 → 5/4) while `del` doubles
(48 → 98) and `loh` rises (18 → 28) — `TRUE` is relabelling and fusing
segments, not resolving new focal events. `test_multi_allelic()` does fire
(it reports finding `13a`, `19g`, `6b`, `11e`, … in some rounds) and reports
`0 multi-allelic CNVs found` in many others; either way the canonical yield
does not improve.

**Recommendation: leave `multi_allelic = FALSE`.** It is what production already
does. The cross-round union (§2) delivers the recovery this parameter was
supposed to, for free and without re-running anything.

## 3b. `min_cells = 10` blocks small samples

`SRX10031191` fails in 19 s with `comparison (==) is possible only for atomic and
list types`. It has ~6 cells after numbat's overlap/SNP filters; `run_numbat()`
does `purrr::keep(subtrees, x$size > min_cells)`, which empties `subtrees`, so
`make_group_bulks()` returns a frame with no `sample` column and `sample` binds
to `base::sample`. Production hit the identical wall on 2026-06-04 — this is not
caused by the rerun. `run_numbat.R` hardcodes `min_cells = 10` at the call site;
lowering it is the only route to a call on that sample, and 6 cells is likely
too few to be worth it.

## 4. `t` — worth a controlled test, but the existing evidence is confounded

`numbat_t: 1e-2` is 1000× numbat's 1e-5 default. Higher `t` means more frequent
HMM state transitions, i.e. a bias toward short segments — which works against
RB's whole-arm drivers, since a fragmented arm splits its evidence and each piece
must clear `min_LLR` alone.

17 samples happen to have run at `t = 1e-5`, so there is an observational
comparison. At **matched `max_iter = 2`**, comparing cross-round stability:

| group | n | mean events, final | mean events, union | gap | zero-event samples |
|---|---|---|---|---|---|
| `t = 1e-2` | 16 | 2.31 | 2.81 | **0.50** | 3 |
| `t = 1e-5` | 17 | 2.88 | 2.94 | **0.06** | 0 |

The `t = 1e-5` runs are nearly stable round to round; the `t = 1e-2` runs lose
half an event per sample between rounds. That is the mechanism §2 documents,
pointing at `t`.

**But this is not a clean comparison and should not be reported as one.** 15 of
the 17 `t = 1e-5` samples are SRR, i.e. the legacy machine's runs — so `t` is
almost perfectly confounded with study, era, and pipeline version. They also ran
`min_LLR = 2` rather than 5, which inflates their call counts independently.

So: a controlled sweep is justified, and this is the highest-value remaining
experiment, but the observational signal above is a reason to run it, not a
result. Note `t` also feeds `detect_clonal_loh()` at `run_numbat.R:259`, so a
sweep moves the LOH prior too — correct behaviour, worth recording.

## 5. Deprioritized: `max_iter` homogenization

The cohort splits 35 samples at 2 rounds / 34 at 4 / 2 at 3, which does confound
cross-sample comparison. But the 2-round samples are almost entirely SRR — 37 of
the 39 rerunnable SRX samples are already at `max_iter = 4`. So this cannot
actually be fixed for the samples where it matters, and the rerun pins the two
stragglers as a side effect. Do not raise it above 4: rounds 5–8 were tested and
bought only wander at ~3.5 h/sample.

## 6. Settled — do not spend compute here

- **`diploid_chroms` pinning.** Tested, negative. Verifiably took effect, did not
  reduce round-to-round wander, weakened two calls, cost 13% runtime.
  (`doc/diploid_baseline_pin_6p.md`)
- **HLA/MHC masking for 6p.** numbat hardcodes an MHC exclusion
  (chr6:28.51–33.48 Mb, 248 genes, all 27 `HLA-*`) inside `filter_genes()`.
  Already cohort-wide. 6p calls are not HLA artifacts. (`doc/hla_masking_6p.md`)
- **`check_convergence = TRUE`.** Tests exact segment equality; never fired in 8
  iterations on either rerun arm.
- **`max_entropy`, `min_cells`.** Already more permissive than default (0.7 vs
  0.5; 10 vs 50).
- **The reference profile.** `/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds`
  (22,901 × 24 retinal cell types) is sound. Beware `data/sridhar_ref.rds` — a
  different, Cone-only file with the same basename that no run has used.

## 7. Reconcile the two `min_LLR` values

Runs use `min_LLR = 5`; `pipeline/scripts/plot_numbat_output.Rmd:103` filters at
50. A 10× discrepancy between calling and plotting should be a stated choice.
Recommendation: keep 5 at the run, and report an LLR 2–5 candidate tier qualified
by the cross-round stability score from §2, rather than lowering the run-time
threshold — which would only add more of the threshold flicker §2 is cleaning up.

## 8. Secondary, untested, low prior

`min_overlap` (0.45 — governs how segments merge into consensus, so it bears
directly on whether an arm survives as one call), `skip_nj = TRUE`, `init_k`,
`min_genes`. Worth looking at only after §4.

---

## How to run it

### Step 1 — the free recovery, all 71 samples

```bash
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript src/diag_rb_scna_recovery.R          # results/rb_scna_*.csv
Rscript src/diag_numbat_param_regimes.R      # results/numbat_param_regimes.csv
```

Both are read-only over `output/numbat_sridhar/`, take under a minute, and touch
no targets store. Already run; outputs are in `results/`.

### Step 2 — `t` pilot, 6 SRX samples, two arms

Pick 6 samples spanning the recovery range — e.g. `SRX10264524` (1 event at
final, 5 in union), `SRX14116944`, `SRX22868104`, `SRX10831283`, `SRX11133590`
(0 events), `SRX10031192`.

```bash
sbatch --array=0-5 --export=ALL,ARM_NAME=ma                       rerun_numbat_arm.sbatch
sbatch --array=0-5 --export=ALL,ARM_NAME=ma_t1e5,T_OVERRIDE=1e-5  rerun_numbat_arm.sbatch
```

`rerun_numbat_arm.sbatch` discovers the 39 rerunnable samples, replays each one's
own logged parameters, and changes only `multi_allelic` (plus `t` in the second
arm), writing to `output/numbat_<ARM_NAME>/`. It never touches
`output/numbat_sridhar/`. To target specific samples rather than the first six
alphabetically, pass the matching array indices — the list is
`ls -d output/numbat_sridhar/SRX*/ | sort`.

Score the arms with the same metric:

```bash
Rscript src/diag_rb_scna_recovery.R output/numbat_ma      _ma
Rscript src/diag_rb_scna_recovery.R output/numbat_ma_t1e5 _ma_t1e5
```

Budget ~3.5 h/sample at `max_iter = 4`, 8 cores, 128 GB.

### Step 3 — full rerun of the 39 SRX samples at the winning arm

> **Not executed — no arm won.** The pilot (§3a) scored 9 union events against
> production's 10, so there is nothing to promote. Do **not** repoint the targets
> pipeline (`numbat_sridhar_dirname` in `pipeline/config.yaml`, or the
> `numbat_rds_files` target) at `output/numbat_ma/`. The command is kept below
> only as the template for a future arm.

```bash
sbatch --export=ALL,ARM_NAME=<arm> rerun_numbat_arm.sbatch   # array 0-38%8
```

Rescore with `src/diag_rb_scna_recovery.R`, compare against the production
**union** — never against production's final round, which is the confound that
made §3's original table look like a win — and only repoint the pipeline once a
new arm demonstrably dominates.

### Step 4 — job-finish ping

```bash
sbatch --dependency=afterany:<JID> --job-name=slack_notify \
  --account=cobrinik_1090 --partition=epyc-64 --time=00:02:00 --mem=1G \
  --output=logs/slack_notify_%j.log \
  --wrap='ST=$(sacct -j <JID> -n -o State | head -1 | xargs); \
          curl -s -X POST -H "Content-type: application/json" \
          --data "{\"text\":\"numbat arm <JID> finished: $ST\"}" "$(cat ~/.slack_webhook)"'
```

All of the above runs outside the targets store, so it does not collide with the
one-`tar_make()`-at-a-time rule. Do not start a `tar_make()` against
`_targets_r431` while step 3 is running.
