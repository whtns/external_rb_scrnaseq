# mosaicMPI NMF gene-program supplement (issue #40)

## Methods
Consensus NMF (cNMF, via mosaicMPI) was run on the 7 paper-retained tumors (SRX10264523, SRX10264526, SRX11133592, SRX11133593, SRX11133594, SRX22868103, SRX10831287).
Raw gene counts per tumor were factorized over a rank sweep, programs were consensus-
aggregated (postprocess), and programs were integrated across tumors into a community
network (single, SCNA-agnostic integration). Recovered programs were compared to the
signatures the main analysis uses: HALLMARK_HYPOXIA, Tirosh S-phase and G2/M markers
(Seurat::cc.genes.updated.2019), and the giotti cell-cycle gene sets. Comparison used
(i) ssGSEA NES of each representative program against those reference sets, aggregated
per community, and (ii) Spearman correlation of each community's per-cell usage
(normalize=False) with the per-cell hypoxia_score / S.Score / G2M.Score.

## Result
Community/communities best matching HALLMARK_HYPOXIA: 1, 4, 5, 13, 14.
Community/communities best matching TIROSH_G2M: 3.
See results/fig_mosaic_nmf_supplement.pdf and results/table_mosaic_nmf_programs.csv.
