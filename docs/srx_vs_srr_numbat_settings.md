# Why the HPC SRX numbat runs differ from the workstation SRR runs

Question: the `SRX*` numbat runs produced on CARC disagree with the earlier `SRR*`
runs produced on the lab workstation. Which settings actually differ?

All statements below are read out of the per-sample `log.txt` parameter block that
numbat writes into `output/numbat_sridhar/<sample>/`, plus `done.txt` (which records
the reference path) and the git history of `pipeline/config.yaml` and
`pipeline/scripts/run_numbat.R`. Nothing here is inferred from the code alone.

## 0. The two sets are the same libraries

`data/metadata.tsv` maps `Run` -> `Experiment` **1:1**. `SRR13884246` and
`SRX10264523` are the same library under a run- vs experiment-level accession; the
rename happened in commit `be5c9fb` ("replace SRR accession IDs with SRX experiment
IDs across R/ and src/"). So the pairs below are directly comparable — no merging of
runs, no extra reads.

The `output/cellranger/<SRX>/outs/` trees are dated **2022-10-14**, i.e. the
alignments were copied over from the workstation rather than re-run on CARC.
Cell Ranger output is *not* a source of difference. (The config now names
`refdata-gex-GRCh38-2024-A`, upgraded from `-2020-A` in commit `be5c9fb`, but no
SRX sample was actually realigned with it.)

`done.txt` gives a clean provenance split:

- every `SRR*` run wrote `/dataVolume/storage/.../sridhar_ref.rds` -> **workstation**
- every `SRX*` run wrote `/project2/cobrinik_1090/.../sridhar_ref.rds` -> **CARC**

Note this means the 2026-03/04/05 `SRR*` runs (numbat 1.3.3 / 1.5.2) were *also*
workstation runs, not HPC ones. "SRR vs SRX" is genuinely "workstation vs HPC".

## 1. The settings that differ

Ordered by how much they can move the calls.

| setting | workstation (SRR) | HPC (SRX) | effect |
|---|---|---|---|
| **numbat version** | 1.2.2 (most), 1.3.0, 1.3.3, one 1.5.2 | **1.5.2** everywhere | different HMM/consensus code, not just parameters |
| **`t`** (HMM transition prob) | `1e-5` on the 2023 runs, `1e-2` from 2024 on | `1e-2` | 1000x. The dominant knob for how many segments get called |
| **`alpha`** | `1e-4` (numbat default — the old call site did not pass `alpha` at all) | **`1e-3` on 31 of 39 samples** | more permissive allele-HMM calls |
| **`min_LLR`** | `2` on the 2023 runs, `5` from 2024 on | `5` | evidence threshold per event |
| **`gamma`** | `50` on part of the 2023 batch, else `20` | `20` | overdispersion prior |
| **`max_entropy`** | `0.8` / `0.5` / `0.7`, varies by sample | `0.7` | how many cells get dropped from clone assignment |
| **`max_iter`** | `2` (one sample `1`) | **`4`** (changed in `be5c9fb`, 2026-06-03) | two extra consensus rounds; segments keep moving after round 2 |
| **`ncores`** | 1–6, varies | 8 | no effect on the result, but affects the NNI search order |
| **cell prefilter** | none | `min_allele_depth = 5`, `min_snps_per_cell = 50` | see below |
| **cell ceiling** | old runs capped at **10 000 cells** | uncapped (`cell_ceiling` is commented out of `run_numbat.R`) | 4 samples were truncated on the workstation |

Unchanged across both: `min_cells = 10`, `max_nni = 100`, `tau = 0.3`,
`skip_nj = TRUE`, `common_diploid = TRUE`, `genome = hg38`, `check_convergence = FALSE`,
`diploid_chroms` unset, and **`multi_allelic = FALSE`** in every production run on
both machines.

`segs_loh` reads as `Given` in the new logs and is absent from the 1.2.2 logs, but
1.2.2 simply did not echo that field — don't read that row as a real change.

## 2. Why the cell counts drop

Commit `162eefd` (2026-03-24) added two prefilters to `pipeline/scripts/run_numbat.R`
that no workstation run ever applied:

- drop allele rows with `n_ref + n_alt < 5` (`min_allele_depth`)
- drop cells with fewer than 50 surviving SNPs (`min_snps_per_cell`)

Every SRX pair loses cells relative to its SRR twin — e.g. SRR13884246 14 131 ->
SRX10264523 12 257, SRR17960480 1 866 -> SRX14116948 1 045. The effect runs the other
way for the four samples the workstation had truncated at 10 000 cells
(SRR13884247/48/49, SRR17960483), which now run at 11–13 k.

Different cell sets change the expression baseline and the initial hclust, so the
clone partition can differ even where every numbat parameter matches.

## 3. The `alpha` trap

`pipeline/config.yaml` currently says `alpha: 1e-4`, with `# alpha: 1e-3` commented
out above it. That is **not** what the production SRX batch ran under:

- `14d58c0` (2026-06-05) set `alpha: 1e-3`
- the main SRX batch ran **2026-06-12** -> `alpha = 1e-3`
- `b4b2dba` (2026-07-13) reverted to `alpha: 1e-4`

So re-running the pipeline from the current config reproduces neither batch. The
eight SRX samples at `alpha = 1e-4` are the ones that ran outside the 06-12 window
(SRX10031191–94 on 06-03, before the change; SRX10831287 and SRX22868103 on 06-15,
via the `rerun_numbat_*.sbatch` scripts, which hardcode their own arguments).

`rerun_numbat_arm.sbatch` already notes this ("config says alpha 1e-4; 31 of 39 SRX
samples ran at alpha 1e-3") and works around it by replaying each sample's parameters
out of its own `log.txt` rather than out of the config.

## 4. Practical reading

The workstation SRR set is not one baseline — it is four years of drifting settings
(t 1e-5 -> 1e-2, min_LLR 2 -> 5, gamma 50 -> 20, numbat 1.2.2 -> 1.5.2) with per-sample
`max_entropy` and a 10 k cell cap. The HPC SRX set is internally consistent except
for `alpha`. Most of the disagreement between a given pair is expected to come from,
in rough order: the `t = 1e-5 -> 1e-2` jump and `min_LLR 2 -> 5` for the 2023-era
samples, the numbat 1.2.2 -> 1.5.2 version gap, `max_iter 2 -> 4`, and the cell-set
change from the new prefilters.

To make a pair genuinely comparable, replay the old parameters from the old
`log.txt` under numbat 1.5.2 — the pattern `rerun_numbat_arm.sbatch` already uses —
rather than comparing the stored outputs directly.

## Appendix: per-pair parameter table

`ws` = workstation (`SRR*`), `HPC` = CARC (`SRX*`).

| pair | run | date | version | t | alpha | gamma | init_k | max_iter | min_LLR | max_entropy | ncores | cells |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SRX10264517 | HPC `SRX10264517` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 6241 |
| | ws `SRR13884240` | 2023-07-31 | 1.2.2 | 1e-05 | 1e-04 | 50 | 5 | 2 | 2 | 0.8 | 6 | 7011 |
| SRX10264518 | HPC `SRX10264518` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 6377 |
| | ws `SRR13884241` | 2023-08-01 | 1.2.2 | 1e-05 | 1e-04 | 50 | 5 | 2 | 2 | 0.8 | 6 | 6847 |
| SRX10264519 | HPC `SRX10264519` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 4085 |
| | ws `SRR13884242` | 2024-03-18 | 1.2.2 | 0.01 | 1e-04 | 20 | 6 | 1 | 5 | 0.5 | 4 | 4197 |
| SRX10264520 | HPC `SRX10264520` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 2109 |
| | ws `SRR13884243` | 2023-03-25 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 2148 |
| SRX10264521 | HPC `SRX10264521` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 7384 |
| | ws `SRR13884244` | 2023-07-31 | 1.2.2 | 1e-05 | 1e-04 | 50 | 5 | 2 | 2 | 0.8 | 6 | 7639 |
| SRX10264522 | HPC `SRX10264522` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 7464 |
| | ws `SRR13884245` | 2026-03-24 | 1.3.3 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 4 | 7464 |
| SRX10264523 | HPC `SRX10264523` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 12257 |
| | ws `SRR13884246` | 2023-09-13 | 1.2.2 | 1e-05 | 1e-04 | 50 | 5 | 2 | 2 | 0.8 | 1 | 14131 |
| SRX10264524 | HPC `SRX10264524` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 11643 |
| | ws `SRR13884247` | 2023-09-18 | 1.2.2 | 1e-05 | 1e-04 | 20 | 6 | 2 | 2 | 0.5 | 4 | 10000 |
| SRX10264525 | HPC `SRX10264525` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 12627 |
| | ws `SRR13884248` | 2023-09-18 | 1.2.2 | 1e-05 | 1e-04 | 20 | 3 | 2 | 2 | 0.5 | 6 | 10000 |
| SRX10264526 | HPC `SRX10264526` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 13194 |
| | ws `SRR13884249` | 2023-09-18 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 10000 |
| SRX11133585 | HPC `SRX11133585` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 1665 |
| | ws `SRR14800543` | 2024-09-23 | 1.2.2 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 3 | 1687 |
| SRX11133586 | HPC `SRX11133586` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 2004 |
| | ws `SRR14800542` | 2026-04-20 | 1.3.3 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.8 | 2 | 2004 |
| SRX11133587 | HPC `SRX11133587` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 2652 |
| | ws `SRR14800541` | 2024-09-23 | 1.2.2 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 3 | 2663 |
| SRX11133588 | HPC `SRX11133588` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 1220 |
| | ws `SRR14800540` | 2024-09-23 | 1.2.2 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 3 | 1244 |
| SRX11133589 | HPC `SRX11133589` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 3184 |
| | ws `SRR14800539` | 2024-09-23 | 1.2.2 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 3 | 3195 |
| SRX11133590 | HPC `SRX11133590` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 1495 |
| | ws `SRR14800538` | 2026-04-20 | 1.3.3 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.8 | 2 | 1493 |
| SRX11133591 | HPC `SRX11133591` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 7013 |
| | ws `SRR14800537` | 2026-03-31 | 1.3.3 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 1 | 7013 |
| SRX11133592 | HPC `SRX11133592` | 2026-06-10 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 10000 |
| | ws `SRR14800536` | 2023-09-18 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 10000 |
| SRX11133593 | HPC `SRX11133593` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 5107 |
| | ws `SRR14800535` | 2023-03-25 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 5481 |
| SRX11133594 | HPC `SRX11133594` | 2026-06-10 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 10000 |
| | ws `SRR14800534` | 2023-09-18 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 10000 |
| SRX14116944 | HPC `SRX14116944` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 7078 |
| | ws `SRR17960484` | 2023-03-25 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 7393 |
| SRX14116945 | HPC `SRX14116945` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 11538 |
| | ws `SRR17960483` | 2023-07-31 | 1.2.2 | 1e-05 | 1e-04 | 50 | 5 | 2 | 2 | 0.8 | 6 | 10000 |
| SRX14116946 | HPC `SRX14116946` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 9442 |
| | ws `SRR17960482` | 2026-05-21 | 1.5.2 | 0.01 | 1e-04 | 20 | 5.0 | 2 | 5 | 0.7 | 1 | 11054 |
| SRX14116947 | HPC `SRX14116947` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 1699 |
| | ws `SRR17960481` | 2023-03-25 | 1.2.2 | 1e-05 | 1e-04 | 20 | 5 | 2 | 2 | 0.5 | 4 | 1833 |
| SRX14116948 | HPC `SRX14116948` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 1045 |
| | ws `SRR17960480` | 2023-08-01 | 1.2.2 | 1e-05 | 1e-04 | 50 | 5 | 2 | 2 | 0.8 | 6 | 1866 |
| SRX22868102 | HPC `SRX22868102` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 5.0 | 4 | 5 | 0.7 | 8 | 9692 |
| | ws `SRR27187902` | 2024-02-03 | 1.3.0 | 0.01 | 1e-04 | 20 | 5 | 2 | 5 | 0.5 | 1 | 10439 |
| SRX22868103 | HPC `SRX22868103` | 2026-06-15 | 1.5.2 | 0.01 | 1e-04 | 20 | 4.0 | 4 | 5 | 0.7 | 8 | 14115 |
| | ws `SRR27187901` | 2024-03-08 | 1.2.2 | 0.01 | 1e-04 | 20 | 4 | 2 | 5 | 0.5 | 1 | 15408 |
| SRX22868104 | HPC `SRX22868104` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 4.0 | 4 | 5 | 0.7 | 8 | 14439 |
| | ws `SRR27187900` | 2024-02-02 | 1.3.0 | 0.01 | 1e-04 | 20 | 4 | 2 | 5 | 0.5 | 1 | 15718 |
| SRX22868105 | HPC `SRX22868105` | 2026-06-12 | 1.5.2 | 0.01 | 0.001 | 20 | 3.0 | 4 | 5 | 0.7 | 8 | 11443 |
| | ws `SRR27187899` | 2024-02-03 | 1.3.0 | 0.01 | 1e-04 | 20 | 3 | 2 | 5 | 0.5 | 1 | 12945 |

---

# Appendix B: the three tumors flagged as thesis-vs-current discrepancies

Three tumors were raised as specific disagreements between the thesis calls and the
current numbat output. All three are reproduced in the stored outputs. ID decoding
from `data/metadata.tsv`:

| shorthand | workstation | HPC | study | sample |
|---|---|---|---|---|
| tumor 524 / "SRR 4247" | `SRR13884247` | `SRX10264524` | wu | RB05_mRNA (GSM5139859) |
| tumor 525 | `SRR13884248` | `SRX10264525` | wu | RB06_mRNA (GSM5139860) |
| tumor 947 / "SRR 481" | `SRR17960481` | `SRX14116947` | field | GSM5883728 |

All three thesis-era runs share one regime — numbat 1.2.2, `t = 1e-5`,
`min_LLR = 2`, `max_entropy = 0.5`, `alpha` unset (1e-4) — and all three HPC runs
share another: numbat 1.5.2, `t = 1e-2`, `alpha = 1e-3`, `min_LLR = 5`,
`max_entropy = 0.7`, `max_iter = 4`.

## Tumor 524 (RB05) — the cleanest test case

The thesis call was **6p gain with 11q loss**. The workstation output on disk
reproduces exactly that:

| round | 6p | 11q |
|---|---|---|
| `segs_consensus_1.tsv` | amp 8.5–40.5 Mb (31.9 Mb), LLR 310 | del 108.2–129.9 Mb (21.6 Mb), LLR 22 |
| `segs_consensus_2.tsv` | amp 8.5–41.1 Mb (32.6 Mb), LLR 321 | del 108.2–129.9 Mb (21.6 Mb), LLR 29 |

One contiguous ~32 Mb 6p block, stable across both rounds, plus a clean 11q deletion.

The HPC run never produces either:

| round | 6p | 11q |
|---|---|---|
| 1 | amp 21.2–33.7 (12.6 Mb) + amp 35.0–36.3 (1.3 Mb) | — (only 11p bits) |
| 2 | amp 0.2–12.9 (12.6 Mb) + amp 17.8–25.5 (7.7 Mb) | `11a` del/bdel 0.2–2.9 Mb — 11**p**, 2.7 Mb |
| 3 | amp 0.2–12.9 (12.6 Mb) + amp 21.2–25.3 (4.1 Mb) | — |
| 4 | none (only a 6q del at 167–170.6 Mb) | `11a` del 0.2–0.6 Mb — 11**p** |

**11q loss is absent from all four HPC rounds.** 6p gain is present in three of four,
but as a *different pair of fragments each round* — which is why
`results/rb_scna_canonical_by_sample.csv` reports it as
`called_final = FALSE, best_stability = 0.25, tier = recovered_marginal` despite
`rounds_called = 1,2,3`. Under the majority-of-rounds policy
(`docs/numbat_round_policy.md`) no individual segment reaches majority, so the arm
drops out of the final call even though the signal is there in three rounds.

Note 11q_loss is not in the canonical event set tracked by the QC tables
(`13q_loss, 16q_loss, 1q_gain, 2p_gain, 6p_gain`), so its disappearance is invisible
in `rb_scna_canonical_by_sample.csv` and only shows up in the raw
`segs_consensus_*.tsv`.

## Tumor 525 (RB06) — the 2p event

Workstation: **2p amp 66 Mb, LLR 222 (round 1) / 65 Mb, LLR 362 (round 2)** — one of
the strongest calls in that sample.

HPC: **no 2p segment in any of the four rounds.** 2p_gain does not appear in
`rb_scna_canonical_by_sample.csv` for `SRX10264525` at all. 6p gain, by contrast,
survives strongly (16 Mb, LLR ~960, all four rounds, `called`).

## Tumor 947 — 9q and the weakened 2p

Workstation rounds 1–2: 1 amp (~130 Mb), **2p amp 85 Mb LLR 58–76**, 6p amp,
**chr9 amp ~100 Mb LLR 34–42**, 13q loh + amp.

HPC rounds 1–4: 1q amp, **2p amp 44 Mb at LLR 6.2**, 6p amp (strong, LLR ~940),
13q del + loh. **No chr9 event in any round.**

So 2p survives but at half the extent and an LLR of 6.2 against a `min_LLR` of 5 —
i.e. it is now a marginal call where it used to be a solid one — and the large chr9
amplification is gone entirely. chr9 is not a canonical tracked event either.

## What this points at

The common signature across all three is that **large contiguous events either
fragment or vanish**: 6p 32 Mb -> 1–13 Mb pieces (524), 11q 22 Mb -> nothing (524),
2p 65 Mb -> nothing (525), chr9 100 Mb -> nothing (947), 2p 85 Mb -> 44 Mb at
threshold (947).

An earlier draft of this section attributed that to the `t = 1e-5 -> 1e-2` change.
**That was wrong, or at least not the leading cause.** Checking the existing rerun
arms on disk settles it for tumor 524.

## Tumor 524 is already answered: it is `multi_allelic`

`SRX10264524` has been rerun twice, and **both reruns reproduce the thesis call at
unchanged production parameters** (`t = 1e-2`, `alpha = 1e-3`, `min_LLR = 5`,
numbat 1.5.2):

| arm | 6p | 11q |
|---|---|---|
| thesis (`SRR13884247`) | one amp 8.5–41.1 Mb (32 Mb), LLR 310 / 321 | del 108.2–129.9 Mb (22 Mb), LLR 22 / 29 |
| production (`multi_allelic=FALSE`) | 13 / 8 / 4 Mb fragments at shifting coords, LLR 8–96; none at all by round 4 | absent in all 4 rounds |
| `numbat_multiallelic` (`TRUE`, no pin) | span restored, stable rounds 1–8, but split: amp 8.5–12.9 (4 Mb, LLR 41) / loh 12.9–21.1 (8 Mb, LLR 24) / **amp 21.2–41.1 (20 Mb, LLR 305)** | **del 108.2–129.9 Mb, LLR 26–32, stable rounds 1–8** |
| `numbat_diploidpin` (`TRUE` + pinned baseline) | same, LLR 292–305 on the 20 Mb block, rounds 1–4 | del 108.2–129.9 Mb, LLR 25–31, rounds 1–4 |

**11q returns at identical coordinates.** 6p returns over the thesis interval but as
a stable three-part amp/loh/amp structure rather than the single 32 Mb amp — a
partial match, worth stating as such.

Both recovering arms share `multi_allelic = TRUE` and nothing else; the multi-allelic
arm has no diploid pin, so the pin is not the cause. `check_convergence` also differs
but only controls early stopping, not the content of any round.

### Why this was not noticed

`rerun_numbat_arm.sbatch` records the `multi_allelic = TRUE` arm as *"TESTED
2026-08-20 AND REJECTED ... zero canonical events gained"*. That verdict was scored on
presence/absence of canonical events against the cross-round union, a metric blind to
both effects at issue here:

- **11q_loss is not in the canonical set at all** — `rb_scna_canonical_by_sample.csv`
  carries only `13q_loss, 16q_loss, 1q_gain, 2p_gain, 6p_gain` — so recovering a
  22 Mb deletion scores as nothing.
- **6p_gain already counted as present** in production's union
  (`tier = recovered_marginal`), so going from shifting LLR-8 fragments to a stable
  LLR-300 block also scores as nothing.

This does **not** overturn the rejection — the over-merging on `SRX10264518` is real —
but the arm was never evaluated against the question the thesis comparison asks.

## The replay job

`rerun_numbat_thesis_replay.sbatch` (drafted 2026-09-13, **not yet submitted**) runs
three arms over the three tumors, `--array=0-8%3`, writing only to
`output/numbat_thesis_replay/<arm>/<sample>/`:

| arm | tasks | what it asks |
|---|---|---|
| `ma` | 0–2 | production params + `multi_allelic=TRUE`. Does the 524 recovery generalise to 525 (2p) and 947 (chr9)? 524 is the positive control. |
| `thesis` | 3–5 | each sample's own workstation params (`t=1e-5`, `alpha=1e-4`, `min_LLR=2`, `max_entropy=0.5`, `max_iter=2`, per-sample `init_k` 6/3/5), prefilters off, numbat 1.5.2. Parameters vs the `multi_allelic` switch. |
| `thesis_122` | 6–8 | identical to `thesis` but against **numbat 1.2.2**. Isolates the version effect. |

### numbat 1.2.2 is installed

Installed 2026-09-13 from the CRAN archive tarball
(`cran.r-project.org/src/contrib/Archive/numbat/numbat_1.2.2.tar.gz`, 4.0 MB,
published 2023-02-14) into **`~/R/numbat-1.2.2`**, a separate library reached only by
prepending it to `R_LIBS` in the `thesis_122` arm. The default library
`~/R/x86_64-pc-linux-gnu-library/4.4` still holds 1.5.2 and is untouched — verified
after install — so the production pipeline is unaffected. The GitHub tag `v1.2.2`
(`kharchenkolab/numbat`, commit `4f43703858`) is the alternative source; the CRAN
tarball was preferred as the artifact the workstation would actually have installed
in 2023.

`pipeline/scripts/run_numbat.R` needs **no modification** to run against 1.2.2: every
argument it passes is in 1.2.2's `run_numbat` signature, and both internals it reaches
for (`numbat:::get_bulk`, `numbat:::detect_clonal_loh`) exist there. 1.2.2's defaults
differ (`min_cells` 50, `init_k` 3, `max_entropy` 0.5, `multi_allelic` TRUE) but the
script passes all of them explicitly. The script asserts the resolved numbat version
per arm and aborts on a mismatch, so a silent fallback cannot turn `thesis_122` into a
duplicate of `thesis`.

### What no arm controls

The cell set. `SRR13884247` and `SRR13884248` were run on the workstation against a
random **10 000-cell subsample** under the old `cell_ceiling` — not reproducible, since
no seed was recorded and `cell_ceiling` is commented out of `run_numbat.R`. Setting
`min_allele_depth=0 min_snps_per_cell=0` in both thesis arms removes the prefilters
added in commit `162eefd`, which is the reproducible part; the subsample is not. So
even `thesis_122` can differ from the archived workstation output on cell membership
alone, and a residual difference there is not evidence of a code or parameter effect.

---

# Appendix C: replay results (job 11959992, completed 2026-09-17)

All nine tasks COMPLETED, exit 0:0, 21 min – 4 h 05 each. Outputs in
`output/numbat_thesis_replay/<arm>/<sample>/`.

## Headline

1. **The discrepancy is parameters, not the numbat version.** The `thesis` arm
   (numbat 1.5.2 + workstation parameters) reproduces the thesis call on all three
   tumors, including the single contiguous 6p block on 524.
2. **numbat 1.2.2 changes essentially nothing at matched parameters.** `thesis` and
   `thesis_122` agree on the events in question, round 1 identical to the LLR.
3. **`multi_allelic=TRUE` alone recovers all three missing events at production
   parameters** — the `ma` arm generalised from 524 to 525 and 947.

## The events

| tumor | event | thesis (archived) | `ma` (prod + multi_allelic) | `thesis` (1.5.2) | `thesis_122` (1.2.2) |
|---|---|---|---|---|---|
| 524 | 6p | amp 8.5–41.1 Mb (32 Mb) LLR 310/321 | split: amp 4 Mb / loh 8 Mb / **amp 21.2–41.1 (20 Mb) LLR 291–305** | **amp 8.5–41.1 (33 Mb) LLR 362/323** | **amp 8.5–41.1 (33 Mb) LLR 362/368** |
| 524 | 11q | del 108.2–129.9 Mb LLR 22/29 | del 108.2–129.9 LLR 26–32 | del 108.2–129.9 LLR 25/31 | del 108.2–129.9 LLR 25/30 |
| 525 | 2p | amp 66 Mb LLR 222 / 65 Mb LLR 362 | **amp 15.6–75.7 (60 Mb) LLR 187–399, all 4 rounds** | amp 9.9–75.7 (66 Mb) LLR 243; 60 Mb LLR 186 | amp 9.9–75.7 (66 Mb) LLR 243; 60 Mb LLR 185 |
| 947 | 2p | amp 85 Mb LLR 58–76 | amp 86 Mb LLR 48; + 17 Mb LLR 32 | amp 0.3–85.6 (85 Mb) LLR 39/78 | amp 0.3–85.6 (85 Mb) LLR 39 |
| 947 | chr9 | amp ~100 Mb LLR 34–42 | **amp 37.7–138.1 (100 Mb) LLR 31–45, all 4 rounds** | amp 37.8–138.1 (100 Mb) LLR 24/43 | amp 37.8–138.1 (100 Mb) LLR 24; r2 64 Mb LLR 12 |

Production called **none** of 524's 11q, 525's 2p, or 947's chr9 in any round.

## Two independent routes back to the thesis calls

- **Thesis parameters** (`t=1e-5`, `alpha=1e-4`, `min_LLR=2`, `max_entropy=0.5`,
  `max_iter=2`, per-sample `init_k`), with `multi_allelic=FALSE` as the thesis runs
  had it. Reproduces the archived workstation output closely — 524 round 2 is
  6p LLR 323 vs the workstation's 321, and 11q LLR 31 vs 29.
- **`multi_allelic=TRUE` at production parameters.** Recovers every missing event,
  but segments 524's 6p as amp/loh/amp over the same interval instead of one block.

So the original `t`-fragmentation hypothesis was right after all for the *contiguous
block*; the `multi_allelic` finding in Appendix B is a second, independent lever. Both
restore 11q, 2p and chr9.

## The version effect is real but small, and not the cause

Jaccard over all non-neutral segments, `thesis` (1.5.2) vs `thesis_122` (1.2.2),
identical parameters:

| sample | round 1 | round 2 |
|---|---|---|
| SRX10264524 | 0.875 (14/15 shared) | 0.350 |
| SRX10264525 | 0.600 | **1.000** |
| SRX14116947 | 0.778 | 0.182 |

Round 1 agrees closely and the tracked events are identical to the LLR. Round 2
diverges on *other* segments (chr1 boundaries, chr15/16/17, and a chr13 `del` vs
`loh` state flip at identical coordinates). So 1.2.2 -> 1.5.2 is **not** a no-op in
general — it just does not explain these three discrepancies.

## The cell set did not matter either

| sample | workstation | `ma` | `thesis` / `thesis_122` |
|---|---|---|---|
| SRX10264524 | 10 000 (random cap) | 11 643 | 11 836 |
| SRX10264525 | 10 000 (random cap) | 12 627 | 12 686 |
| SRX14116947 | 1 833 | 1 699 | 1 829 |

`ma` matches production exactly, confirming it differs from production only by
`multi_allelic`. The thesis arms ran on ~12 k cells against the workstation's random
10 k subsample for 524/525 and still reproduced the calls — so the unreproducible
subsample flagged in Appendix B turned out not to matter for these events. For 947
the cell set is reproduced to within 4 cells.

---

# Appendix D: cohort-wide t=1e-5 sweep (job 12134445, completed 2026-09-20)

36 tasks: 34 COMPLETED, 1 FAILED (`SRX10031191`, the known `min_cells` wall),
1 OUT_OF_MEMORY (`SRX22868103`, 128 GB at 4h16 — `t=1e-5` costs more memory because
segments are longer). Runtimes 42 min – 4h43. Output `output/numbat_t1e5/`.

Only `t` differs from each sample's production run; every other parameter was
replayed from that sample's own `log.txt`, and `multi_allelic` was held at the
production value `FALSE`.

## 1. `t` controls fragmentation, as predicted — measured on 34 matched samples

| arm | `t` | segs/round | median seg | focal <10 Mb | total segs |
|---|---|---|---|---|---|
| production | 1e-2 | 22.8 | 10.2 Mb | 50% | 3004 |
| sweep | 1e-5 | 15.3 | 32.4 Mb | 25% | 1980 |

`t=1e-2` calls **49% more segments at a third the median length**, half of them focal.

## 2. Canonical events, 32 matched samples, identically rescored

Production was rescored with the current script (`_prodnow`) so both sides are
comparable; it reproduced the August numbers exactly (172 final / 209 union on the
full 71-sample set), confirming the scoring is stable.

| arm | FINAL | UNION | cross-round gap |
|---|---|---|---|
| production | 82 | 113 | **+31 (38%)** |
| `t=1e-5` | **104** | 106 | **+2 (2%)** |

| event | prod final | t1e5 final | prod union | t1e5 union |
|---|---|---|---|---|
| 13q_loss | 17 | 17 | 23 | 17 |
| 16q_loss | 17 | **21** | 22 | 21 |
| 1q_gain | 22 | **27** | 26 | 27 |
| 2p_gain | 12 | **20** | 20 | 21 |
| 6p_gain | 14 | **19** | 22 | 20 |

**`t=1e-5`'s final round (104) beats production's final round (82) by 27%, and nearly
matches production's cross-round union (113) without needing the union at all.**
Production is slightly ahead on union total (113 vs 106), driven by 13q and 6p.

## 3. Segment quality — the axis presence-counting misses

Canonical-event segments, 32 matched samples:

| arm | n segs | median extent | median stability | called in all rounds | median LLR |
|---|---|---|---|---|---|
| production | 451 | 18.7 Mb | 0.25 | 10% | 69.4 |
| `t=1e-5` | 271 | **46.6 Mb** | **0.50** | **27%** | **217.0** |

Fewer segments, 2.5× the extent, 2× the stability, 3× the likelihood ratio.

## 4. Per-sample, it almost never loses

| | samples |
|---|---|
| `t=1e-5` better | 16 |
| tie | 15 |
| `t=1e-5` worse | **1** (`SRX10264522`, loses 13q_loss) |

## 5. Recommendation

**Set `numbat_t: 1e-5`** in `pipeline/config.yaml` — numbat's own default. It is
better or equal on 31 of 32 samples, raises final-round canonical events 27%, and
removes the cross-round instability that round selection exists to work around.

Consequences to plan for:

- **Round selection becomes largely redundant.** Its value was capturing a 22–38%
  cross-round gap; at `t=1e-5` that gap is 2%. Do not remove the machinery, but the
  selected round and the final round will mostly coincide.
- **`SRX22868103` needs ~256 GB**, not 128.
- Production's union still catches ~7 events `t=1e-5` misses (13q and 6p mainly).
  A `t=1e-5` + union combination is untested.

## 6. What this is not

There is still no orthogonal SCNA truth for these tumors, so this is a
stringency-matched comparison, not an accuracy measurement. Longer segments with
higher LLR are what you expect from *correct* whole-arm calls, but also from
over-merging; the two are not distinguished here. The five-event canonical frame is
also unchanged, and 11q — the event that started this whole investigation — is still
outside it.

## 7. CORRECTION (2026-09-21): three samples drop out entirely at t=1e-5

The sweep's "34 COMPLETED" is misleading. Three samples exit 0 and write `done.txt`
while producing no usable numbat object, because every candidate CNV is filtered away:

| sample | numbat's stop message | production rounds |
|---|---|---|
| `SRX11133590` | No CNV remains after filtering by **LLR in pseudobulks** (suggests lowering `min_LLR`) | 4 |
| `SRX10031192` | No CNV remains after filtering by **entropy in single cells** (suggests raising `max_entropy`) | 2 |
| `SRX11133586` | No CNV remains after filtering by **entropy in single cells** | 2 |

Coherent mechanism: longer segments make per-cell posteriors more diffuse, so cells
fail the `max_entropy = 0.7` gate. This is a real cost of `t=1e-5` that the
canonical-event comparison in §2 does not show, because these samples are absent from
the t1e5 scoring and were therefore excluded from **both** sides of the matched
comparison. The 82 -> 104 headline is unaffected, but the honest per-sample tally over
the 37 attempted SRX samples is:

- 16 better, 15 tie, 1 worse (`SRX10264522`, loses 13q_loss)
- **3 produce nothing usable**
- 1 unrunnable at any `t` (`SRX10031191`, the `min_cells` wall)
- 1 needed 240 GB (`SRX22868103`, OOM at 128 GB)

`t=1e-5` remains the recommendation, but it is not free: it trades three whole samples
for better segmentation on the rest. If those three matter, they need their own
`max_entropy` / `min_LLR` treatment rather than a cohort-wide setting.

## 8. Round selection becomes unnecessary at t=1e-5 (measured)

`src/select_numbat_round.R output/numbat_t1e5 _t1e5`:

```
samples where selection beats final: 1   unchanged: 33   worse: 0
selected round still misses this many union events: 0
mean non-neutral segments — selected 15.2, final 14.9
```

At `t=1e-2` round selection was worth +34 canonical events (172 -> 206). At `t=1e-5`
it is worth **one sample**. Use the final round and drop selection — which also
removes the selecting-on-the-outcome bias, since the per-sample maximum over rounds
is not an unbiased estimate.
