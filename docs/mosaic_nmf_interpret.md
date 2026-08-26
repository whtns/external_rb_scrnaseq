# mosaicMPI NMF -- unbiased interpretation (issue #40 add-on)

## Purpose
Companion to the supplement: instead of projecting communities onto our 4 chosen
hypoxia/cell-cycle references, characterise every community data-driven -- against
the full 50-set MSigDB HALLMARK library (P1), from its own top marker genes (P2),
against our Seurat clusters (P3), and spatially on each tumor's UMAP (P4).

Pan-tumor communities (present in >=5/7 tumors): 1, 2, 3, 4.

## Top HALLMARK match per community
| community | top HALLMARK | NES |
|---|---|---|
| 1 | OXIDATIVE_PHOSPHORYLATION | 0.256 |
| 2 | E2F_TARGETS | 0.498 |
| 3 | G2M_CHECKPOINT | 0.398 |
| 4 | HYPOXIA | 0.176 |
| 5 | INTERFERON_ALPHA_RESPONSE | 0.479 |
| 6 | NOTCH_SIGNALING | 0.295 |
| 7 | P53_PATHWAY | 0.350 |
| 8 | INTERFERON_ALPHA_RESPONSE | 0.382 |
| 9 | E2F_TARGETS | 0.498 |
| 10 | MYC_TARGETS_V2 | 0.271 |
| 11 | MYC_TARGETS_V2 | 0.253 |
| 12 | MYC_TARGETS_V1 | 0.220 |
| 13 | NOTCH_SIGNALING | 0.235 |
| 14 | MYC_TARGETS_V2 | 0.321 |
| 15 | E2F_TARGETS | 0.369 |
| 16 | TGF_BETA_SIGNALING | 0.268 |

See results/fig_mosaic_nmf_interpret.pdf and results/table_mosaic_nmf_community_genes.csv.
