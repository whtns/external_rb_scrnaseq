# RB scRNA-seq Targets Pipeline Overview

Pipeline for single-cell analysis of retinoblastoma tumors with somatic copy-number
alterations (SCNAs: 1q+, 2p+, 6p+, 16q−). Built with the `{targets}` framework.

---

```mermaid
flowchart TD

    %% ── CONFIGURATION ───────────────────────────────────────────────────────
    subgraph CFG["⚙️ Configuration"]
        direction TB
        IS[interesting_samples\n20 SRR samples]
        LCC[large_clone_comparisons\nclone pairs per sample]
        LCS[large_clone_simplifications\narm → segment map]
        LFE[large_filter_expressions\nCNV filter logic]
        CD[cluster_dictionary]
        BO[branch_dictionary\nCO[cluster_orders]]
        SM[subtype_markers\nLiu & Radvanyi 2022]
        MP[mps\n3CA MetaPrograms]
        IG[interesting_genes\n1q / 2p / 6p / 16q / S1 / S2]
        RB[rb_scna_samples\nper-SCNA sample lists]
    end

    %% ── INPUT FILES ─────────────────────────────────────────────────────────
    subgraph IN["📂 Input Files"]
        direction TB
        NRF[numbat_rds_files\nper interesting sample]
        SEUS[seus\nSeurat objects per sample]
        ORIG[original_seus\nunfiltered Seurat objects]
        DBF[debranched_seu_files\npre-computed file paths]
    end

    %% ── STUDY METADATA ──────────────────────────────────────────────────────
    subgraph META["📋 Study Metadata & QC"]
        direction TB
        SCS[study_cell_stats]
        FS04[fig_s04\nStudy cell stats]
        TS01[table_s01\nUMI / genes / mito %]
        TS02[table_s02\nSample metadata]
        TS03[table_s03\nFiltered cluster tally]
        FS02T[fig_s02 / table_s04\nTCGA SCNA frequency]
        FS03[fig_s03\nGISTIC plot]
    end

    %% ── SEURAT PROCESSING ───────────────────────────────────────────────────
    subgraph SP["🔬 Seurat Processing"]
        direction TB
        USEU[unfiltered_seus\nCNV-labelled cells]
        FSEU[filtered_seus\ncluster + cell QC]
        FINSEU[final_seus\nset_final_seus]
        RSEU[regressed_seus\nCC-regressed]
        HYPSEU[hypoxia_seus\nhypoxia-scored]
        LOWH[seus_low_hypoxia\nhypoxia ≤ 0.5]
        HIGHH[seus_high_hypoxia\nhypoxia > 0.5]
    end

    %% ── DEBRANCHED SCNA SUBSETS ─────────────────────────────────────────────
    subgraph DB["🌿 Debranched SCNA Subsets"]
        direction TB
        DBSEU[debranched_seus\nall 23 branches]
        DB1Q[debranched_seus_1q\n5 samples]
        DB2P[debranched_seus_2p\n4 samples]
        DB6P[debranched_seus_6p\n4 samples]
        DB16Q[debranched_seus_16q\n3 samples]
        SCNA[scna_seus\nunion of all SCNA sets]
        HYP1Q[hypoxia_seus_1q/2p/6p/16q\nhypoxia-filtered subsets]
    end

    %% ── CNV / NUMBAT ANALYSIS ───────────────────────────────────────────────
    subgraph CNV["🧬 CNV / Numbat Analysis"]
        direction TB
        LNPDF[large_numbat_pdfs\nraw Numbat plot files]
        FLPF[filtered_large_plot_files\nCNV-filtered plots]
        FS03A[fig_s03a\nunfiltered CNV heatmaps]
        FS05[fig_s05\nsmoothed expression]
        FS13[fig_s13\nfiltered CNV heatmaps]
        FS06A[fig_s06a\nkaryograms]
    end

    %% ── DIFFERENTIAL EXPRESSION ─────────────────────────────────────────────
    subgraph DE["📊 Differential Expression"]
        direction TB
        ADE[all_diffex_clones\nall genes, per branch]
        CDE[cis_diffex_clones\nin-segment genes]
        TDE[trans_diffex_clones\nout-of-segment genes]
        ADEC[*_diffex_clones_for_each_cluster\ncluster-stratified versions]
        TADE[table_all_diffex_clones\nExcel output]
    end

    %% ── INTEGRATION ─────────────────────────────────────────────────────────
    subgraph INT["🔗 Integration by SCNA"]
        direction TB
        I1Q[integrated_seus_1q\n4 samples]
        I16Q[integrated_seus_16q\n3 samples]
        I2P[integrated_seus_2p\n2 samples]
        I6P[integrated_seus_6p\n2 samples]
        ILOWH[integrated_seu_*_low_hypoxia\nhypoxia-stratified integrations]
    end

    %% ── CORRESPONDING STATE ANALYSIS ────────────────────────────────────────
    subgraph CS["🔄 Corresponding State Analysis"]
        direction TB
        CSEU[corresponding_seus\nmatched sample pairs]
        CSD[corresponding_states_dictionary]
        CSDE[corresponding_clusters_diffex\nclone-matched diffex]
        CSEN[corresponding_clusters_enrichments]
        SD2P[states_dictionary_2p/6p\nSCNA-specific state maps]
        CSDE2P[corresponding_clusters_diffex_2p/6p]
    end

    %% ── ENRICHMENT & ONCOPRINT ──────────────────────────────────────────────
    subgraph EN["🔍 Enrichment & Oncoprint"]
        direction TB
        OENR[oncoprint_enrich_clones_gobp\nGO-BP enrichment]
        OENH[oncoprint_enrich_clones_hallmark\nHallmark enrichment]
        OENRH[oncoprint_enrich_clusters_hallmark\ncluster-level Hallmark]
        UOI[unfiltered_oncoprint_input_by_scna]
        OI[oncoprint_input_by_scna\nfiltered oncoprint data]
        OP[oncoprint_plots\nmain oncoprint figures]
    end

    %% ── VISUALIZATION ───────────────────────────────────────────────────────
    subgraph VIZ["🎨 Visualization"]
        direction TB
        CT[clustrees\nper sample]
        CTC[clustree_compilation PDF]
        AHC[annotated_heatmap_collages]
        CC[collage_compilation PDF]
        DCT[debranched_clone_trees]
        SVL[subtype_violins]
        FH[factor_heatmaps_1q/2p/6p/16q]
        PPP[patchwork_phase_plot]
        CLD[clone_pearls_1q/2p/6p/16q]
        CLSC[cluster_scorecard]
    end

    %% ── MANUSCRIPT FIGURES ──────────────────────────────────────────────────
    subgraph FIG["📄 Manuscript Figures"]
        direction LR
        F02[fig_02\nIntegrated 1q all samples]
        F03[fig_03\n16q- integration]
        F04[fig_04\n2p+ integration]
        F05[fig_05\n6p+ integration]
        F09[fig_09\n2p+ cluster DE]
        F10[fig_10\n6p+ cluster DE]
        FS07[fig_s07\n1q sample-specific]
        FS11[fig_s11\n2p supplemental]
        FS12[fig_s12\n16q sample-specific]
        FS0408[fig_s04_08\n2p integration]
        F07[fig_07a/b\n1q cluster DE ± integration]
        F08[fig_08a/b\n16q cluster DE ± integration]
        FS08[fig_s08/s09/s10/s20\nCluster-level DE supp. figs]
    end

    %% ── FINAL OUTPUTS ───────────────────────────────────────────────────────
    subgraph OUT["✅ Final Outputs"]
        FAT[figures_and_tables\naggregated list]
        PN[pipeline_notification\nemail on completion]
    end

    %% ── EDGES ───────────────────────────────────────────────────────────────

    IN  --> SP
    CFG --> SP
    SP  --> DB
    DB  --> CNV
    DB  --> DE
    CFG --> DE

    SP      --> INT
    DB      --> INT
    INT     --> CS
    CSEU    --> CS
    CS      --> EN
    DE      --> EN
    CFG     --> EN

    EN  --> FIG
    CNV --> FIG
    INT --> FIG
    CS  --> FIG
    DB  --> FIG

    DB  --> VIZ
    DE  --> VIZ
    EN  --> VIZ

    META --> OUT
    FIG  --> OUT
    VIZ  --> OUT
    FAT  --> PN
```

---

## Pipeline Stages Summary

| Stage | Key Targets | Description |
|-------|------------|-------------|
| **Configuration** | `interesting_samples`, `large_clone_*`, `cluster_dictionary` | Constants and per-sample clone/segment mappings |
| **Input Files** | `numbat_rds_files`, `seus`, `debranched_seu_files` | Raw Numbat and Seurat RDS file paths |
| **Study Metadata** | `study_cell_stats`, `fig_s02–s04`, `table_s01–s03` | QC metrics, TCGA SCNA frequency, GISTIC |
| **Seurat Processing** | `filtered_seus` → `final_seus` → `hypoxia_seus` | Filtering, QC, CC regression, hypoxia scoring |
| **Debranched SCNA Subsets** | `debranched_seus_1q/2p/6p/16q`, `scna_seus` | Clone-branch–specific Seurat objects per SCNA type |
| **CNV / Numbat** | `large_numbat_pdfs`, `fig_s05`, `fig_s13` | Numbat CNV heatmaps, expression roll plots, karyograms |
| **Differential Expression** | `all/cis/trans_diffex_clones` (± per cluster) | Clone-vs-clone DE genome-wide, cis, and trans |
| **Integration** | `integrated_seus_1q/2p/6p/16q` | Cross-sample integration by SCNA type |
| **Corresponding States** | `corresponding_clusters_diffex_2p/6p` | Matched clone-state DE across integrated samples |
| **Enrichment & Oncoprint** | `oncoprint_input_by_scna`, `oncoprint_plots` | GO-BP / Hallmark enrichment, oncoprint figures |
| **Visualization** | `clustrees`, `annotated_heatmap_collages`, `factor_heatmaps_*` | Clustrees, heatmap collages, clone trees, violin plots |
| **Figures** | `fig_02–10`, `fig_s07–s25` | All manuscript and supplementary figures |
| **Final Output** | `figures_and_tables`, `pipeline_notification` | Aggregated output list + email notification |

---

## Key Dependencies (simplified)

```
interesting_samples ──► debranched_seu_files ──► debranched_seus ──► [1q/2p/6p/16q subsets]
                                                                             │
                     ┌───────────────────────────────────────────────────────┤
                     ▼                     ▼                    ▼            ▼
              CNV heatmaps          Diffex clones         Integration    Collages
              (fig_s03a/s05/s13)    (all/cis/trans)       (per SCNA)     (fig_s07/s12...)
                     │                     │                    │
                     └──────────►  Enrichment / Oncoprints ◄───┘
                                          │
                                    Manuscript figures
                                    (fig_02–10, fig_s*)
                                          │
                                   figures_and_tables
                                          │
                                  pipeline_notification
```

## Notes on Broken / WIP Targets

| Target | Issue |
|--------|-------|
| `collages_1q` and hypoxia variants | Reference `hypoxia_1q` / `hypoxia_1q_low` / `hypoxia_1q_high` which are **not defined**. Commented out with TODO. Use `hypoxia_seus_1q` once fixed. |
| `heatmap_collages_hypoxia_low/high` | Were **identical copies** of `heatmap_collages_hypoxia`. Commented out with TODO to use `seus_low_hypoxia` / `seus_high_hypoxia`. |
| `fig_s11`, `fig_s23` | Call `not_sure_what_this_does()` — **placeholder function name** needs to be replaced. |
| `fig_s04_09` | Output path was `"results_fig_s04_09.pdf"` (missing `/`). Fixed to `"results/fig_s04_09.pdf"`. |
