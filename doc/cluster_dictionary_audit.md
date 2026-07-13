# Cluster-dictionary filtering audit

_Generated 2026-06-30 from `results/cluster_dictionary/cluster_dictionary_discrepancies.tsv` + `cell_metadata` cell counts._

Quantifies cells removed by the auto cluster_dictionary (clusters flagged
`remove=1`, dropped wholesale as `cluster_remove` in `annotate_filter_reason`).

## Samples ranked by % cells removed (dictionary only)

| sample | total cells | clusters removed | cells removed | % removed |
|---|---:|---:|---:|---:|
| SRX10831280 | 2991 | 2 | 2477 | 82.82 |
| SRX10831281 | 529 | 1 | 346 | 65.41 |
| SRX11133588 | 1580 | 2 | 767 | 48.54 |
| SRX14116948 | 1866 | 3 | 877 | 47 |
| SRX11133585 | 2083 | 1 | 855 | 41.05 |
| SRX14116945 | 12452 | 1 | 5005 | 40.19 |
| SRX10264524 | 11836 | 1 | 3621 | 30.59 |
| SRX14116944 | 7393 | 1 | 2072 | 28.03 |
| SRX10264525 | 12686 | 1 | 3494 | 27.54 |
| SRX10264521 | 7596 | 1 | 1731 | 22.79 |
| SRX14116947 | 1833 | 2 | 416 | 22.7 |
| SRX10264520 | 2140 | 2 | 460 | 21.5 |
| SRX10264522 | 7638 | 1 | 1312 | 17.18 |
| SRX10264523 | 12335 | 3 | 2039 | 16.53 |
| SRX22868103 | 15408 | 2 | 2357 | 15.3 |
| SRX22868102 | 10439 | 1 | 1517 | 14.53 |
| SRX22868104 | 15718 | 2 | 1680 | 10.69 |
| SRX10264518 | 6847 | 3 | 701 | 10.24 |
| SRX10264517 | 7011 | 2 | 533 | 7.6 |
| SRX11133593 | 5481 | 1 | 409 | 7.46 |
| SRX14116946 | 9490 | 1 | 692 | 7.29 |
| SRX11133592 | 10798 | 1 | 577 | 5.34 |
| SRX10264526 | 13400 | 2 | 618 | 4.61 |
| SRX11133589 | 3257 | 2 | 65 | 2 |
| SRX11133594 | 12195 | 1 | 220 | 1.8 |
| SRX10031191 | 1449 | 0 | 0 | 0 |
| SRX10031192 | 8967 | 0 | 0 | 0 |
| SRX10031194 | 2550 | 0 | 0 | 0 |
| SRX10264519 | 4197 | 0 | 0 | 0 |
| SRX10831282 | 5276 | 0 | 0 | 0 |
| SRX10831283 | 5214 | 0 | 0 | 0 |
| SRX10831286 | 990 | 0 | 0 | 0 |
| SRX10831287 | 9153 | 0 | 0 | 0 |
| SRX11133586 | 4031 | 0 | 0 | 0 |
| SRX11133587 | 2785 | 0 | 0 | 0 |
| SRX11133590 | 3152 | 0 | 0 | 0 |
| SRX11133591 | 7354 | 0 | 0 | 0 |
| SRX22868105 | 11522 | 0 | 0 | 0 |

## Program totals (clusters / cells removed)

| program | clusters | cells |
|---|---:|---:|
| low_qual | 14 | 16288 |
| rod | 11 | 2894 |
| MALAT1 | 7 | 11047 |
| cone | 5 | 4274 |
| APOE | 3 | 338 |

## Low-confidence removals (3) — review these

rod/cone flagged on only the 2-gene rule minimum, or incidental `top1_gene`
removals (e.g. APOE). These are the likeliest false-positive deletions.

| sample | cluster | abbrev | rule | n_matched | cells | % | top_genes |
|---|---:|---|---|---:|---:|---:|---|
| SRX11133594 | 5 | APOE | top1_gene | 1 | 220 | 1.8 | APOE, B2M, VIM, PDC, SAT1, HLA-C, HLA-B, FOS, CEBPD, YBX3 |
| SRX10264518 | 9 | APOE | top1_gene | 1 | 80 | 1.17 | APOE, VIM, CLU, GPX3, PTN, LGALS1, TF, IFITM3, FRZB, B2M |
| SRX11133589 | 4 | APOE | top1_gene | 1 | 38 | 1.17 | APOE, B2M, SPP1, FTL, HLA-DRA, APOC1, CD74, VIM, SAT1, SRGN |
