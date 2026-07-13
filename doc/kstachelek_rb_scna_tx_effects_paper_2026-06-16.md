

**PNAS Template for Main Manuscript**

This PNAS template for the Main Manuscript may be used to organize your main text source file. The template is intended to provide a clearly organized PDF to facilitate the review process. Further information is available in our [Author Center](https://www.pnas.org/authors/submitting-your-manuscript#article-types).

**Using the template**

Paste the appropriate text from your manuscript into the relevant section of the template. You may maintain the template formatting, or reapply styles after pasting your text into the template.

Figures should be placed on separate pages with legends set immediately below each figure. Table titles should be set immediately above each table. Please do not use field codes for figures and tables to ensure the numbering is correct after the PDF conversion.

References cited in the main text should be included in a separate reference list at the end of this file. Examples of the PNAS citation style are included below.

**Notes about submission**

The following items are required on the title page: **Title, Author Line, Author Affiliations, Corresponding Author information**. Each of these items **must** be provided on the title page for us to proceed with processing your paper.

You are not required to adhere to the section order outlined below. For example, you may combine your Results and Discussion or use alternate section headings. Materials and Methods should be included after the Results and Discussion in most cases. If your paper does not include the standard section headings, please provide a brief explanation in the **Comments to Editorial Staff** field of the submission form.

You may include subheadings within the standard headings listed below. You may also include line numbers.

**Submitting your main manuscript**

Delete this first page, and then save your completed main text file as a PDF for submission, following instructions available [here](http://www.pnas.org/site/authors/procedures.xhtml#preparation).

6 pages, but PNAS allows articles up to a maximum of 12 pages. A standard 6-page article is approximately 4,000 words, 50 references, and 4 medium-size graphical elements (i.e., figures and tables).

*Updated February 2022*

See emails:

Re: hypoxia score 10/10/25, 3:24 pm and reply at 3;47 pm and KS paper notes 10/9/25 11:24 am

![](data:image/png;base64...)

**Main Manuscript for**

**Chromosome arm-level somatic copy number alterations are associated with altered cell cycle distributions and novel cell states in retinoblastoma tumors**

Kevin Stachelek1,2\* and David Cobrinik1,3,4,5,\*

1 The Vision Center, Department of Surgery, and Saban Research Institute, Children's Hospital Los Angeles, Los Angeles, CA, USA.

2 Cancer Biology and Genomics Program, Keck School of Medicine, University of Southern California, Los Angeles, CA, USA.

3 Department of Cancer Biology, Keck School of Medicine, University of Southern California, Los Angeles, CA, USA.

4 Norris Comprehensive Cancer Center, University of Southern California, Los Angeles, CA, USA.

5 USC Roski Eye Institute, Department of Ophthalmology, Keck School of Medicine, University of Southern California, Los Angeles, CA, USA.

\*Correspondance: David Cobrinik and Kevin Stachelek

**Email:** stachele@usc.edu, \* dcobrinik@chla.usc.edu

**Author Contributions:** Paste the author contributions here.

**Competing Interest Statement:** The authors declare no competing interests.

**Classification:** Paste the major and minor classification here. Dual classifications are permitted, but cannot be within the same major classification.

**Keywords:** Somatic copy number alterations (SCNAs), retinoblastoma, tumorigenesis, single-cell RNA-sequencing, (ONE MORE?)

**This PDF file includes:**

Main Text

Figures 1 to X

Tables 1 to X

# Abstract

Most retinoblastomas initiate in response to biallelic *RB1* loss and progress by acquiring secondary mutations, resulting in tumors composed of genomically heterogeneous subclones. Chromosome arm-level somatic copy number alterations (SCNAs) 1q+, 2p+, 6p+ and 16q- are the most common secondary changes, yet their associated cell behaviors are unclear. Here, we present a generalizable approach to define cell state changes associated with recurrent SCNAs that i) uses Numbat analysis of retinoblastoma single cell RNA-sequencing (scRNA-seq) datasets to delineate SCNA-defined subclones, ii) separately defines hypoxic and non-hypoxic cell state distributions of tumor subclones that differ by an SCNA of interest, and iii) identifies candidate SCNA driver genes with cell-state-specific gene expression changes. This approach revealed that subclones with 16q- were depleted in a post-mitotic G0 state and enriched in G1 and S phase states with increased hypoxia markers and decreased expression of metallothionein genes located on chromosome 16q. Subclones with 1q+ in addition to 16q- had increased proportions of cells in G1- S-, and G2-like states and a novel G2/M-like state with increased expression of S/G2 marker genes located on 1q. A subclone with 2p gain including the *MYCN* oncogene (but not one lacking MYCN) had diminished contribution to a differentiated-cone-like state and increased contribution to a G1-like state with upregulated MYCN targets. Partial 6p gains [Gains of proximal and distal 6p regions] were associated with distinct cell state changes and candidate cancer progression drivers ~~depending on the gained 6p region~~. These findings show that recurrent RB SCNAs are associated with cell state redistributions, acquisition of novel states, and cell-state-specific gene expression changes consistent with retinoblastoma progression.

**Significance Statement**

Most cancers have chromosome arm-level gains and losses called somatic copy number alterations (SCNAs). Highly recurrent SCNAs are thought to contribute to tumorigenesis, yet SCNA-associated changes in cancer cell behavior have been unclear. This study presents an approach SCAnState (SCNA-associated Cell Analysis of State Transitionse) to identify cell-state changes that are associated with recurrent SCNAs in the childhood eye cancer, retinoblastoma. The study reveals that acquisition of recurrent retinoblastoma SCNAs is accompanied i) by a redistribution of tumor cells between the cell states that are present in the ancestral non-SCNA-containing tumor cells and ii) by the appearance of new cell-states, with each SCNA associating with different cell state changes. Defining SCNA-associated cell state changes may enable new therapies that counteract SCNA oncogenic effects.

# Introduction

Cancers develop in response to successive genetic and epigenetic changes in a susceptible cell-of-origin (Greaves & Maley 2012). Cancer-related genetic changes include single nucleotide variants (SNVs) and focal copy number alterations, usually affecting known cancer genes with well-defined effects, and larger whole-chromosome and chromosome-arm-level somatic copy number alterations (SCNAs), usually encompassing hundreds of genes and having poorly defined effects (Shlien & Malkin 2009; Ben-David & Amon 2020). These large SCNAs can appear in premalignancies (Mitchell et al. 2018; Körber et al. 2019; Douville et al. 2021, Dimaras 2008), during emergence of a malignant clone (Erickson et al. 2022; Girish et al. 2023), or during cancer progression, after tumors emerge (Sen et al. 2022; Martins et al. 2022). In aggregate, they occur in nearly 90% of solid tumors and in 50% of hematopoietic cancers, with an average of 10 arm-level SCNAs per malignancy (Shukla et al. 2020).

Each type of cancer harbors a unique pattern of recurrent SCNAs, suggesting that cancer-specific recurrent SCNAs harbor oncogenes (for gains) or tumor suppressor genes (for losses) that contribute to tumorigenesis (Taylor et al. 2018). Indeed, the frequency of chromosome arm gains and losses across many cancers correlates with the number of oncogenes and tumor suppressor genes on each arm (Davoli et al. 2013; Taylor et al. 2018). As such, specific SCNA effects may be attributed to the combined gain or loss of many genes rather than a single gene on affected chromosomes (Xue et al. 2012; Bonney et al. 2015). The combined gain or loss of genes is thought to enhance cancer cell fitness (Laughney et al. 2015; Shih et al. 2023; Mennie et al. 2026), yet it has been challenging to discern how they do so.

One approach to investigating SCNA effects is to compare overall (*i.e.*, *bulk*) gene expression in histologically similar tumors with versus without an SCNA of interest. However, such comparisons may be confounded by intra-tumor heterogeneity for the SCNA, by the variable presence of additional genomic changes, by host-specific genetic background effects, and by tumor microenvironment features, like hypoxia, which are extrinsic to the cancer’s genetic determinants (Bhattacharya et al. 2020). Moreover, tumor cells exist in many distinct gene expression states related to cell differentiation, cell cycle phase, and cell stress (Tirosh & Suva 2024), which may or may not be affected by the SCNA under study. If most tumor cells in a sample exist in states that are unaffected by an SCNA, the copy-number-related gene expression changes that drive tumorigenesis may be difficult to detect (Figure 1A).

Single cell analysis of tumors with subclonal SCNAs may overcome these challenges. Simultaneous capture of each cell’s DNA and RNA content through single cell DNA + RNA sequencing or by inference of DNA copy number from single-cell RNA sequencing (scRNA-seq) (Vandereyken et al. 2023; Sdeor et al. 2024) enables exploration of cancer cell states both within and across the distinct subclones in individual tumors. However, while several single-cell sequencing studies have sought to relate SCNA-associated gene expression to tumorigenesis (Tirosh & Suva 2024), this approach has been limited by the low throughput of single cell DNA + RNA sequencing and the somewhat imprecise copy number assignments of initial scRNA-seq-based methods. To overcome these limitations, methods that integrate haplotype information with allele and gene expression signals, such as Numbat or XClone, enable more precise identification of SCNA-defined tumor subclones (Gao et al. 2022; Huang et al. 2024; Song et al. 2025)). However, to our knowledge, such methods have not been used to delineate SCNA-associated cell state changes.

The current study uses Numbat analysis of scRNA-seq datasets to define the cell state changes associated with highly recurrent cancer-progression-related SCNAs in the childhood retinal cancer, retinoblastoma. Retinoblastoma is an apt context in which to study SCNA-associated cell state changes, as most retinoblastomas have the same cancer-initiating insult and cell-of-origin (biallelic *RB1* loss in cone photoreceptor precursors), carry one or more of the highly recurrent ‘RB SCNAs’ (1q+, 2p+, 6p+, 16q-), and have low secondary mutation burdens (an average of four SNVs in protein-coding genes (Kooi et al. 2016a, 2016b; Gröbner et al. 2018; Stachelek et al. 2023; Cobrinik 2024)), which limits the confounding effects of co-occuring genomic changes. Moreover, although RB SCNAs are occasionally seen in putative premalignant regions (Dimaras 2008, Eagle), most appear to be acquired after retinoblastomas emerge and are subclonal (Kooi et al. 2016), which enables intratumor subclone-specific gene expression comparisons. Candidate RB SCNA driver genes have been proposed based on their differential expression in bulk RNA-sequencing of retinoblastomas with and without RB SCNAs of interest (Kooi 2016 PLOS1, Supplementary Fig. **S1**). However, as for other such comparisons, it is unclear if the altered expression of such genes drives the expansion of SCNA-containing subclones or reflects secondarily altered cell state distributions (Figure 1A). Additionally, most retinoblastomas have hypoxic and non-hypoxic regions (Sudhakar et al. 2013; Zhang et al. 2025; Wan et al. 2025), whose effects on SCNA-driven cell state changes is currently poorly understood (ref). Thus, understanding RB-SCNA-associated cell state changes in retinoblastoma may also reveal SCNA driver effects in the other malignancies where these SCNAs frequently appear (Fig. **S1**,**S2,** Table **S1**).

To elucidate RB SCNA-associated transcriptomic cell state changes, we developed an scRNA-seq-based approach in which i) SCNA-defined tumor subclones are identified, ii) tumor cells are partitioned according to a hypoxia gene expression index, iii) hypoxic and non-hypoxic cell states are separately annotated according to cell cycle and cell differentiation features, iv) cell state distributions are determined for each tumor subclone, and v) genes that are differentially expressed between subclones in the same or highly similar cell states are identified (Figure 1b). This approach reveals that acquisition of RB SCNAs is associated with a redistribution of the tumor cells among the already-existing cell states, with the appearance of novel cell states, and with cell state-specific differential expression of candidate SCNA drivers.

# Results

### Retinoblastoma scRNA-seq datasets with subclonal acquisition of RB SCNAs

To evaluate RB SCNA-associated cell state changes, we retrieved and examined thirty-one publicly available retinoblastoma scRNA-seq datasets from five studies (Collin et al. 2021; Yang et al. 2021; Wu et al. 2022; Field et al. 2022; Liu et al. 2024) (Table S2). Seven samples were excluded because they failed Numbat analysis (likely due to [high subclone complexity or?] sample quality) and an additional x# were excluded due to their high proportion of cells with low numbers of unique reads or genes detected or their high percentages of mitochondrial reads (Figure S3, Table S2) in accord with scRNA-seq best practices (Luecken & Theis 2019)). Exceptions were made for two [or three?] samples ([SRX10031194?] SRX11133588 and SRX11133585) whose high mitochondrial reads were attributed to contaminating normally mitochondria-rich photoreceptors.

For the remaining 21 datasets, we performed Numbat analysis (Gao et al. 2022) and Louvain clustering to identify and remove clusters comprised of non-tumor cells based on their lack of SCNAs and their assigned non-tumor cell identities.on high-quality cells, defined as having >103 genes detected and <5% mitochondrial reads. Numbat processes scRNA-seq data to detect all SCNAs in a sample, assigns a probability that each cell harbors each of the identified SCNAs, and groups cells with highest probability of carrying the same SCNAs into subclones (Figure 1b) (Extended Data 1, Table S3). To ensure the robustness of downstream analyses, we excluded cells in which the highest probability of any SCNA assignment was <50%. This revealed clusters composed exclusively of diploid cells while marker gene analyses assigned each cluster’s cell type (Figure S4). Diploid clusters were removed from the dataset except those annotated as cones, which in most cases were considered to be early retinoblastoma cells derived from the cone precursor cell of origin. However, clusters composed of diploid cone-like cells were not removed in two [or three?] samples that had extensive normal retina contamination([SRX10031194, ?] SRX11133588 and SRX11133585), in which wild type cones are likely present. (Table S4). Apart from these highly retina-contaminated samples, an average of 2.4% of high quality cells were excluded per sample.

A custom script generated a) segmentation plots with cells grouped according to the identified SCNA subclones, b) plots of each cell’s probability of having the SCNAs identified in each sample, c) phylogenetic relationships between the subclones, d) karyograms displaying gained and lost chromosome segments, and e) smoothed gene expression profiles of cells clustered based on copy number inferred from gene expression alone, which is similar to the output of inferCNV and can be used to cross-check Numbat results (Fig. 2, Fig. S5-7). Seven samples had complex aneuploidy that did not yield clear phylogenetic relationships (Supplementary Fig. S5), while five had interpretable phylogenies yet lacked subclone pairs that differed by an RB SCNA of interest (Supplementary Fig. S6). These samples were removed from further analysis, although their diploid tumor cells were included in the integrated diploid RB subclone (see below).

Nine samples contained one or more pairs of subclones that differ by one or two RB SCNAs, including four with acquired 1q+, four with acquired 2p+, four with acquired 6p+, and three with acquired 16q- (Fig. clone\_tree\_collage). Eleven samples with insufficient cell numbers in relevant subclones, widespread chromosomal instability or multiple co-acquired SCNAs were excluded (Fig. clone\_tree\_collage). The reproducibility of the inferred SCNAs was demonstrated by the similar SCNA phylogenies and segmentation plots for each of three pairs of scRNA-seq replicates (Supp Fig. clone\_tree\_collage\_of\_merged\_replicates). Although four tumors had sequential or concurent 16q- and 1q+ as initial SCNAs (Fig. 1b), we observed no significantly preferred order of SCNA acquisition based on what test, consistent with prior evaluations of bulk tumors (Kooi et al. 2015, 2016a).

For samples with acquired RB SCNAs, we performed Louvain clustering over a range of resolutions and examined whether different subclones were differentially distributed among the different clusters. In a pilot analysis, tumor SRX11133594 was inferred to have three SCNA subclones (diploid, 16q-, and 16q-, 1q+) that were unevenly distributed among the tumor cell clusters at all resolutions examined (Fig. 2A,B, S#). See Supplementary Information for explanations for choosing optimal resolutions for main figures (*Fig S7, S11, S15*). Then, we visualized relationships between clusters in UMAP plots, which segregated cells in part according to their Louvain clusters and their G1, S, and G2/M assignments determined as described [or – using standard S and G2/M markers (Tirosh et al. 2016)) (Fig. 2C. However, most clusters were composed of cells from more than one cell cycle phase, which hindered interpretation of the cell state redistributions.

To better define cell cycle status, cells were plotted in cell cycle space using standard cell cycle marker genes (Tirosh, Stuart 2019. This enabled ordering of clusters according to average cell cycle position, as marked with centroids that also indicate the proportion of cells in each cluster(Fig. 2D, E).

Comparing the subclone distributions in each cluster revealed that the 16q- subclone had an increased proportion of cells in a post-mitotic state whereas the 16q-.1q+ subclone contributed to a new cluster with an exceptionally high G2/M score (Fig. 2F). For most clusters, the subclone distributions differed from the tumor as a whole ~~based on a binomial test~~. This suggested that acquisition of 16q- and 1q+ was associated with altered distribution of tumor cells among the pre-existing cell states and with the emergence of a new G2/M-like state.

*Let’s hold off on the following, which may go into more detail than needed at this stage, in favor of addressing hypoxia first. Can keep the text here to see if some points can be migrated to the subsequent sections. See Stratification of retinoblasotma cell states based on hypoxia below.*

While the analysis in Fig. 2 displays the different cell cycle distributions of the three SRX11133594 tumor subclones, defining the effects of individual SCNAs requires binary comparison of tumor cell subclones with and without a single SCNA of interest. To enable binary comparisons between subclones, after the initial running of Numbat, we filtered the cells so that only the initial clone with the SCNA of interest and its immediate antecedent clone without the SCNA were retained, and then re-ran Louvain clustering. For example, to evaluate transcriptome changes related to 1q gain in SRX11133594, the SRX11133594 diploid tumor cells were filtered so only the 16q- and 16q-, 1q+ clonesremained, followed by SCNA assignment, and clustering.

To gain insight into the distinct features of each cluster, we plotted heatmaps displaying each cluster’s subclone compositions, S and G2/M scores, and marker gene expression (Fig. 3A, showing only the 16q- and 16q-,1q+ clones for simplicity). Examination of many such plots revealed that the precise marker genes of clusters with similar cell cycle positions varied, that clusters of genomically highly complex tumors with many SCNAs differend from those with few SCNAs, and that retinoblastomas consistently included (*to next p*)

varied for different clusters in different tumors, we distinguished G0 and/or G1 states that differed according to expression of markers of cone photoreceptor differentiation such as *GNB3* (Ritchey et al. 2010) and *GUK1* (Kallman et al. 2020). We also distinguished two S phase clusters in most samples: a major S cluster often marked by high expression of S phase marker genes *HMGB1 (*high mobility group box 1*), RRM2* (ribonucleotide reductase), and various histone genes, and a smaller population distinguished by high expression of S and M-related chromatin regulators *22T* and *ATAD2* and having a higher S phase score than the more abundant S phase cells, which we labeled S+. S/G2 clusters have additional G2-related gene expression including TOP2A, CENPF, and UBE2C. We also frequently identified a post-mitotic (pm) cell cluster characterized by expression of G2 signature with low S score and expression of residual G2 genes including *PTTG1* and *ARL6IP1*. We also identify

transcriptomic states characterized by markers of heat shock (e.g., HSPA1A or HSPB1 (O’Flanagan et al. 2019)), and hypoxia (BNIP3, ENO1, GAS5), with cell cycle positions often spanning g1/s/g2 or centered on g1, respectively (Fig. 3). Thus, cell cycle variation was a major contribution to distinct cell states along with stress states distinguished by hypoxia and heat shock markers. These observations prompted us to integrate retinoblastoma scRNA-seq datasets showing subclonal heterogeneity for the same SCNA to reveal SCNA-related cell state changes common to multiple tumors and to separately explore SCNA-related cell state changes in hypoxic and non-hypoxic cell populations.

**Stratification of retinoblastoma cell states based on hypoxia-related gene expression.**

In the SRX11133594 pilot analysis, most Louvain clusters centered on a specific cell cycle position, yet some clusters spanned multiple positions in cell cycle space and were marked by hypoxia-related genes (*BNIP3*, *ENO1*, *GAS5*) and/or heat shock-related genes (*HSPA1A* or *HSPB1* (O’Flanagan et al. 2019)) (Fig. XX), implying that the clusters were more strongly defined by cell stress rather than by cell cycle phase. All clusters marked by heat shock genes also highly expressed hypoxia-related genes, whereas the converse was not the case, implying that the heat shock response primarily occurred in cells expressing a hypoxia response. …. Notably, some cell stress clusters that spanned several cell cycle phases divided into cell cycle-specific clusters at higher resolution (Fig. S#), which confounded attempts to define SCNA-associated cell state changes.

As hypoxia profoundly affects cell behavior (ref) and potentially alters SCNA effects, we segregated cells according to their expression of a hypoxia-related gene expression index tailored to retinoblastoma scRNA-seq datasets, RB Hypoxia index genes were selected based on their presence in previously developed hypoxia indices (or gene ontologies?) such as *BNIP3*, *ENO1*, and *GAS5 or* other genes marking the same clusters*??*, including …, …, (Figure?). Inspection of multiple retinoblastoma scRNA-seq datasets revealed that all? clusters with high expression of hypoxia-related genes had low expression of mitochondrial H-strand genes (*MT-CO2*, *MT-ATP6*, *MT-CO3*, *MT-ND3*, and *MT-CYB*) (Figure 3a?). The loss of mitochondrial gene expression is consistent with a hypoxic state (ref), although to our knowledge, the selective loss of H strand (but not L strand?) expression is novel. These and other considerations led to our x# - gene hypoxia index incorporating both high expression of hypoxia-related genes and low expression of mitochondrial H strand genes (Supplementary Fig #), which was used to stratify cells into hypoxic and non-hypoxic states.

Addtion of the hypoxia stratification step completed the scRNA-seq based SCNA and cell state profiling workflow: 1) select high-quality scRNA-seq datasets, 2) remove low-quality cells, 3) SCNA profile and remove cells with low confidence SCNA assignments, 4) produce and annotate Louvain clusters and remove clusters composed of diploid non-tumor cells, 5) stratify cells based on hypoxia-related gene expression, 6) separately Louvain cluster and annotate hypoxic and non-hypoxic cell states at several clustering resolutions, 7) define SCNA subclone distributions among the different hypoxic and non-hypoxic cell states at different clustering resolution (Fig. 3b).

*It seems best to treat 16q- and 1q+ together and in that order given all relevant samples have both and they’re aquired in that order. E.g.,*

*1) 16q vs diploid individual tumors*

 *16q vs diploid (combine the individual tumors)*

 *16q vs diploid (combine the individual tumors AND even more diploids from other tumors)*

*2) )16q-,1q+ vs 16q-, individual tumors*

 *16q-,1q+ vs 16q-, combine the individual tumors*

 *16q-,1q+ vs diploid individual,*

 *16q-,1q+ vs diploid combined including SRX10264526 (where there’s no intermediate 16q- wubclone)*

### Tumor subclones with 16q- are depleted in post-mitotic and enriched in hypoxia cell states

We identified ten samples with 16q- including a lone initiating 16q- in three tumors (SRX11133594, SRX11133593, and SRX11133592). We examined three samples with 16q- acquired over a diploid background and integrated tumor clones with 16q- and antecedent diploid clones. We identified cluster marker genes indicative of cell cycle state and characteristic stress states (hypoxia or hsp) (Fig. 5A). In the non-hypoxic group, G1 states were distinguishable by expression of photoreceptor marker genes including *GNB3* and so were labeled G1 PR 1 and 2 whereas cell cycle expression better characterized additional cell states including G1/S with expression of *TYMS*, S with expression of *PCLAF*, G2 with expression of *CENPF* and a post-mitotic state with expression of *ARL6IP1*.

In the hypoxic group, … Cell stress states included a distinct heat shock state with expression of *HSPA1B* and two hypoxia states, both with high expression of *BNIP3* and one with *SOX11* (De Bolòs et al. 2024)and *CCNB2*, the latter consistent with hypoxic cells in a G2/M phase. We then characterized the distribution of clones with and without 16q- across states. The largest significant differences were a reduction in the percentage of cells in a post-mitotic (pm) state in clones with 16q- (6.4%) compared to antecedent diploid clones (11.7%, p<.001) and enrichment of a hypoxia state (12.3% versus 7.4%, p<.001) (Fig. 5C-E; Table S4b). Other transcriptomic states did not show consistent divergence in clone proportion with or without 16q-.

To identify specific effects of 16q- we compared differential expression between clones with and without 16q- within cell states. We considered differentially expressed genes in cis to 16q- SCNA segments and all other differentially expressed genes in trans. We observed differential expression between clones with and without 16q- in cis in G1, G1/S, S, and pm cell states of metallothionein genes (Fig. S10), located in a gene cluster at 16q13, whereas we detected no differential expression between clones in the G2 state. As with comparisons within cell states, comparisons across all clusters in cis showed dramatic variation in metallothionein factors *MT1X and MT1F* (Fig. S10).

Variations in 16q- state were not an artifact of clustering resolution. We also evaluated clone distributions for each tumor separately to assess whether SCNA-related cell state distributions observed in the integrated dataset were independently observed in individual samples (Fig. S12, Fig. S13). In two of three samples, pm state was significantly depleted in the 16q- clone (5.2-7.9%) versus the antecedent clone (7.4-15.7%). Likewise, hypoxia was enriched (8.6-17.2%) versus antecedent (5.8-12.8%) (Fig. S12, Table S4bB). Thus, we observed similar patterns of enrichment of the 16q- clone according to cell cycle state in analysis of each sample.

### Tumor subclones with 1q+ are enriched in S/G2 phases and occupy a novel G2/M state.

We identified 12 tumors with chromosome 1q+, among which one tumor gained 1q together with 16q- and three gained 1q after 16q- (Fig. S24). The remaining eight samples with 1q+ had insufficient 1q diploid cells with which to compare 1q+ cells. We examined 1q+ and 1q+,16q- associated transcriptomic changes in an integrated analysis of these four samples (Fig. 4). Integration was performed solely on the subset of cells in each tumor populating the least evolved 1q+ clone and its antecedent 1q diploid clone.

We performed Louvain clustering of the integrated dataset and plotted the expression of cluster marker genes in a heatmap annotated with tumor, cluster, and SCNA subclone identifiers, along with an empirically-derived G2M.Score and S.Score (Fig. 4A). We next plotted cells colored by Louvain cluster on an x-y axis according to expression of established cell-cycle phase marker genes (Tirosh et al. 2016) and established average cell cycle position by plotting the cluster centroids separately plotting clusters with higher expression of hypoxia markers BNIP3 and GAS5. At the optimal clustering resolution, this revealed seven clusters with well-defined cell cycle positions centered at G0/G1, G1/S, S/G2, G2/M, post-mitotic phases, and three defined by markers of hypoxia, centered either at G0/G1, S, or G2/M (the latter of which also had high-level heat shock protein gene expression. (Fig. 4B, C).

We also assessed whether SCNA-related cell state distributions in the integrated dataset were observed in individual samples (Fig. S5, Fig. S6). In each sample, state G1 MT was significantly enriched in the 1q+ clone (3-14.2%) versus the antecedent clone (1-7.4%). Likewise, S/G2 was enriched (10-15.7%) versus antecedent (3-5.1%) and the G2/M state was more dramatically enriched with 1q+ (9.7-14.6%) versus (0.1-0.9%) (Fig. S5, Table S4b). Thus, we observed similar patterns of enrichment of the 1q+ clone according to cell cycle state in analysis of each sample (Fig. S5, S6).

We distinguished three G0 or G1 states according to expression of markers of photoreceptor differentiation and mitochondrial expression. State G1 DIFF had the lowest S phase score with expression of cone markers including *GNB3* (Ritchey et al. 2010)and *DCT* (Kallman et al. 2020) and the pan-photoreceptor marker RCVRN; state G1 MT showed elevated mitochondrial gene expression; and state G1 SUB2 showed elevated expression of retinoblastoma subtype 2 marker genes *(Liu et al. 2021)* *(*Fig. 4A, Fig. S25). We distinguished two S phase clusters: a major S cluster with high expression of S phase marker genes *CCNE2* which peaks at the G1-S phase transition (Caldon & Musgrove 2010) and *RRM2*, and a smaller population with expression of S and M-related chromatin regulators *SMC4* (cohesin) and *CENPU* which we labeled S+(Fig. 4A)*.* We identified an S/G2 cluster with expression of the S and G2-related gene expression including *TOP2A, CENPF and UBE2C*, while a cluster designatedG2-M had expression of G2 markers without expression of S phase*.* Cells labeled PM (post-mitotic) are characterized by a decreased G2 score and residual expression of some G2 genes including *ARL6IP1*. We also identify transcriptomic states characterized by markers of cell stress including heat shock (hsp) markers such as *HSPA1A*, and hypoxia markers including *BNIP3*,with cell cycle score centroids overlapping G1 or spanning G1/S/G2, respectively (Fig. 4B-C). Notably, the hypoxia cluster had elevated expression of hsp genes and the hsp cluster had elevated expression of hypoxia genes, while both had reduced mitochondrial RNA expression, as did the S+ cluster, suggesting that cells in these clusters experience similar stress yet in different cell cycle phases.

We next characterized the distribution of tumor clones with and without 1q+ across cell states. By measuring the proportion of cells in clones contributing to each transcriptomic state we found striking differences between clones with or without 1q+ which correlated with cell cycle progression.

Notably, a G2/M state was populated almost entirely by 1q+ clones (Fig. 4D). Indeed, 12.2% of cells in the 1q+ clones were assigned to the G2/M cluster, whereas <1% of cells in the antecedent clones were assigned to this cluster (p<.001, binomial test) (Fig. 4E, *integrated plot*). There was a less dramatic imbalance in the S/G2 states (~14.8% of 1q+ versus ~4.5% of antecedent, p<.001) and in the G1 MT state (10.4% v. 5.5%, p<.001). In contrast, a larger proportion of antecedent clones occupied the G1 DIFF state with expression of differentiated photoreceptor marker genes (15.0% of 1q+ versus 28.0% of antecedent) and occupied the G1 SUB2 state (12.3% of 1q+ versus 19.2% of antecedent).

We observed a lower proportion of cells in hypoxia states within clones with 1q+ (11.3%) versus antecedent clones (24.6%, p<0.001)) whereas the proportion of cells in a heat shock state was similarly increased (2.2% v. 1.7%, p<0.05\* for heat shock; 1.6% v. 0.8%, P<.001 for S+) (Fig. 4E, Table S4b). The most abundant G1 state in clones with 1q+ showed elevated mitochondrial gene expression while cells in hypoxic, heat-shock, and S+ states showed notably decreased expression of mitochondrial genes.

### Candidate 1q+ driver genes identified by differential expression are enriched for proliferation-related genes and gene sets

To identify specific effects of 1q+, we compared gene expression between clones with and without 1q+ in each transcriptomic cluster, among genes located on 1q (in *cis*) or outside of 1q (in *trans*). In the four tumors with acquisition of 16q- and 1q+, we observed differential expression in *cis* in G1, S, and S/G2 cell states (Fig. S8A-B). Overexpression ranged from 20.5 to 21.1 (1.4 to 2.1-fold) consistent with the 1.5-fold increase in gene dosage for cells with three rather than two chromosome 1q copies. Surprisingly, there was little overlap between overexpressed 1q genes detected in different G0/G1 clusters and less overlap still between overexpressed S cluster genes and overexpressed G1 cluster genes~~and no overlap between such genes~~. Differential expression in S/G2 is limited to *CEP350, CENPL and ZNF670,* suggesting that the increase in 1q+ proportions may relate 1q+ overexpression in preceding cell cycle phases.

Through analysis of global gene expression changes (in *trans*), the G1 SUB2 cluster showed enrichment of TNFA signaling via NFKB with differential expression of pathway members *ATF3 and IER5* (Fig. S8B). The G1 DIFF cluster likewise showed enrichment of TNFA SIGNALING VIA NFKB, driven by *ATF3, and IER5 in cis.* Interestingly, the G1 MT cluster showed enrichment of G2M CHECKPOINT related genes, including *CENPF,* and *CKS1B* in cis as well as *HCN3* . Notably, *HCN3* cooperates with *MYC* to promote mouse hepatocellular carcinoma and is frequently upregulated in human hepatocellular carcinoma where it is also co-increased with *MYCN* (Zhang et al. 2023). S state cells had enrichment of gene sets including APOPTOSIS and HYPOXIA, with elevated expression of differentially expressed genes including *ATF3, LMNA,* (Fig. S8, Fig. S9).

Differential expression between clones across all clusters yielded significantly differentially expressed genes in cis with dramatic variation in G2 proliferation markers including *CENPF, NUF2, NEK2* (Fig. S8, Fig. S9).), likely due to the enrichment of 1q+ cells in advanced cell cycle states (Fig. 4).

To distinguish G1 states, we performed differential expression between clusters, regardless of subclone composition. Comparing G1 MT and G1 DIFF states (enriched or depleted for 1q+ cells, respectively) showed enrichment (in G1 MT) of proliferation-related hallmark gene sets including *G2M checkpoints* (Fig. S8, Fig. S9) consistent with higher S phase score of the G1 MT state*.* Comparing G1 SUB2 versus G1 DIFF and G1 MT, we also observed enrichment of msigdb oncogenic signature gene sets upregulated upon E2F1 activation (Fig. S9). Thus, the 1q+ enriched G1 MT population has properties similar to G1 DIFF yet has higher S phase score suggesting further cell cycle progression.

### Tumor subclones with 2p+ are enriched in a G1 cell state with activation of G2M checkpoint and E2F targets ontologies

We identified 9 tumors affected by gain of all or part of 2p. One of these (SRX10264523) gained all of chromosome 2 and one (SRX10264524) gained a portion of 2p that lacked MYCN. Among tumors with any 2p+ gain, the SCNA was preceded or accompanied by 1q+ in 6 tumors and by 6p+ in 3 tumors. One tumor showed 2p+ preceded by 1q+, 16q- (SRX10264526, and two tumors showed 2p+ preceded by 6p+,11q- or exclusively by 6p+ (SRX10264524 and SRX10264525). 2p+ did not arise as a unique initiating event in any case (Fig. S14). To identify effects of 2p+ we examined one sample with 2p+ acquired over a 6p+ background and one sample with combined 2p+,10q+ acquired over 1q+,6p+, then integrated 2p+ clones and antecedent clones. In the integrated data set, we observed multiple G1 states, two with low S phase scores (g1a1 and g1a2) and two with higher S phase scores (g1b and g1c). Among these, g1a2 had notable expression of cone markers including *PDC* whereas g1b had notable mitochondrial gene expression*.* Cell cycle states S/G2 and G2/M were characterized by expression of *RRM2* and *TOP2A*, respectively, with a second S/G2 cluster having notably high mitochondrial expression (Fig. S14). A post-mitotic state had high expression of *ARL6IP1.* We observed two cell stress states with expression of hypoxia markers including *ZFAS1* as well as heat shock markers *HSPA1B and HSPB1* (Fig. S14A-B*).*

We next characterized the distribution of tumor clones with and without 2p+ across cell states and found striking differences G1 cell states. In the integrated dataset, G1 state A2 was enriched dramatically in antecedent clones (comprising 7.4% of non-2p+ cells) and depleted in the 2p+ clone (comprising 0.4% of 2p+ cells) (Table S4B, Fig. S15). Similar results were obtained when the same clustering was applied to the individual SRX10264524 and SRX10264525 samples (Fig 4.9E, *arrows*). We also observed SCNA-related cell state distributions in individual samples (Fig. S15, S16B).

To interrogate variation between the four G0/G1 transcriptomic states, we performed differential expression between clusters with diverging clone makeup. Differential expression in trans coinciding with 2p+ showed activation of hallmark gene set MYC TARGETS and downregulation of photoreceptor genes including *PDC* and activation of *TFF1*, a marker of an aggressive retinoblastoma subtype (Fig. S17). We identified candidate driver genes in cis including *POLE4* and *PPM1G* (Fig. S18). Interestingly, *MCYN* expression did not significantly vary in comparisons between subclones, possibly due to low mean expression across all cells.

### Select proposed candidate SCNA driver genes are identified by single cell analysis

We identified genes that were differentially expressed between SCNA-affected and antecedent subclones in all samples and compared this list with putative retinoblastoma SCNA driver genes (Fig. 6). We focused on affected genes in cis (i.e., in the gained or lost segment) that were previously noted in studies of candidate SCNA driver genes (Kooi et al. 2016a). We found 14 previously proposed candidate genes including *KIF14, MCL1, MDM4, ZNF281* in 1q and *CAP2, CCND3, DEK, E2F3, HLA-A, KIF13A, NUP153, PIM1, SOX4, and TDP2* in 6p (Table S8). Overlap of candidate driver genes identified with gene-dosage effects included 1q: *ZNF281*, and 6p: *CCND3, DEK, E2F3, KIF13A, NUP153, PIM1, SOX4, TDP2, and 16q: CYLD (Kooi et al. 2016a).* We also identified several cases where single cell sequencing data supports proposed candidate driver genes that were not verified by prior gene-dosage measurements including *KIF14, MCL1, MDM4* all within 1q segments.

# Discussion

Retinoblastoma tumor progression proceeds by acquisition of secondary genomic changes which contribute to more aggressive and treatment-resistant forms. Secondary mutations apart from *RB1* and *MYCN* in the form ofnucleotide substitution, insertion, or deletion mutations are rare (<0.1 mutations/Mb) (Gröbner et al. 2018). Somatic copy number alterations (SCNAs) likely play a greater role in tumor progression based on their greater frequency in bulk tumor analyses. Relative to other cancers, however, retinoblastoma genomes are remarkably stable.

In the course of the study, we identified a hypoxia-related gene expression signature that enables separate evaluation of SCNA-associated cell state changes in hypoxic and non-hypoxic cell populations

Retinoblastomas display a burden of arm-level SCNAs equivalent to the second lowest of the 33 cancer types present in TCGA (Kooi et al. 2016b). A subset of SCNAs in retinoblastoma (RB SCNAs) are however, highly recurrent including 1q+, 2p+, 6p+, and 16q- (Kooi et al. 2016a).

1q+ is one of the most frequent SCNAs in hepatocellular carcinoma, occurring in 50-70% of cases. This amplification includes genes like MCL1, IL6R, and CHD1L that promote cell survival, inflammation, and genomic instability. 1q+ typically appears early in hepatocarcinogenesis and is associated with more aggressive disease, increased recurrence risk, and poorer survival outcomes (Jacobs & Norton 2021). 2p+ frequently occurs in uterine carcinosarcoma and often overlaps MYCN. The presence of 2p+ often indicates worse clinical outcomes and may represent a potential therapeutic target (Zhao et al. 2016). 16q- is a characteristic feature of low-grade breast neoplasia, including atypical ductal hyperplasia and low-grade ductal carcinoma in situ, suggesting it's an early event in breast carcinogenesis. This deletion targets several tumor suppressor genes including CDH1, CBFB, and CTCF (Harbers et al. 2021).

Retinoblastoma is therefore an appropriate setting to investigate the role of somatic copy number alterations in the process of tumor progression in cancer, generally.

Subclonal heterogeneity for RB SCNAs has been observed in a large fraction of tumors and presents an opportunity to compare tumor subclones with otherwise identical genetic backgrounds. However, most studies of SCNA effects in retinoblastoma have relied on bulk tumor sequencing, which necessitates rigorous control of confounding experimental variables as well as assembly of large tumor cohorts. This work brings together for the first time all published retinoblastoma tumor scRNA-seq data to define associations between recurrent somatic copy number alterations and transcriptomic cell states. Previous analyses of these data identify tumor heterogeneity in cell states defined by known cell type-specific markers, expression cluster markers, and pathway analysis of cluster markers (Field et al. 2022). One study also makes use of methylation, RNA expression, and bulk copy number information to define retinoblastoma subtypes and to sketch subtype heterogeneity within a tumor scRNA-seq sample (Liu et al. 2021). These investigations attempt to define functional drivers of Rb progression and in some cases identify correlations between dedifferentiation, proliferation-related expression, and aneuploidy but do not report analyses of phenotypic effects of specific SCNAs. In contrast, we examined samples with subclonal heterogeneity for each RB SCNA for evidence of associations between SCNA incidence and cell state. We determined whether each of the RB SCNAs enable new and distinct transcriptomic states.

Our approach uses inference of SCNA clones from droplet scRNA-seq data with sufficient power for clonal comparison of clones that have acquired recurrent RB SCNAs. In this work we evaluated the use of Numbat to identify transcriptomic effects of somatic copy number alterations (SCNAs) in retinoblastoma scRNA-seq samples and successfully inferred recurrent ‘RB SCNAs’ (1q+, 2p+, 6p+, and 16q-), in high quality scRNA-seq datasets. Our analyses also revealed retinoblastoma tumor cell states. Most cell clusters were restricted to specific positions in cell cycle space with expression of known photoreceptor and neuronal cell-types in G1, or markers of S and G2 cell cycle phases as well as markers of cell stress including hypoxia and heat shock. Discrimination of cell states enables comparisons of SCNA-associated changes in retinoblastoma tumors. We identified SCNAs that are present in distinct cell clusters and therefore enable new cell states and SCNAs that are enriched or depleted in individual cell states relative to their frequency in the whole tumor. As a caveat, it is also possible that associations between acquired SCNAs and cell states relate to unobserved point mutations or epigenetic changes.

To identify retinoblastoma cell states we sought programs of expression that are reflected in the most common expression programs of intratumor heterogeneity (Gavish et al. 2023). Cancer is characterized by uncontrolled cell growth, with tumors containing actively dividing cancer cells. scRNA-seq provides a static snapshot, typically capturing only a fraction of cancer cells during their cell cycle. The proportion of cycling cells varies widely among tumors, from less than 1% in slow-growing tumors to a majority in rapidly proliferating tumors and cell line models. Cycling cells can be categorized into distinct cell cycle phases, based on S phase and G2/M phase scores. Non-cycling cells may be in a temporary G1 state or in a longer-term quiescent or senescent state (G0). Retaining cell cycle phase information is crucial for understanding tumor characteristics and cell state distributions. Cancer cells also vary by patterns of cellular stress-related programs, including responses to hypoxia, heat shock, and unfolded proteins.

To mitigate the impact of inter-tumor heterogeneity, integration methods designed to remove technical batch effects were employed for four samples with subclonal 1q+, three subclonal 16q-, two with subclonal 2p+. Integration methods may retain some patient-specific signals while removing others, complicating the interpretation of remaining inter-tumor differences. To mitigate distortions due to integration of diverse tumors sharing a common RB SCNA, we chose to integrate exclusively tumor clones with a specific acquired SCNA and antecedent clones. To the extent possible, this method limited confounding tumor variation due to unobserved mutations and additional co-acquired SCNAs. It is also possible that associations between acquired SCNAs and cell states relate to unobserved variations in point mutations. We also assume maximum parsimony in our designation of SCNA evolutionary history when determining the order of SCNA events in cancer progression and in our identification of candidate driver genes.

By examining integrated scRNAseq data from tumors with isolated 1q+, we found a novel G2/M cell state composed entirely of 1q+ tumor subclones. We also investigated whether SCNAs are enriched or depleted in individual cell states relative to their frequency in a tumor as a whole and identified enrichment of 1q+ clones in cell states corresponding to post-mitotic and S/G2 cell cycle phases (Fig. 4). Differential expression between tumor subclones with and without 1q+ in distinct cell states revealed upregulation of hallmark gene sets including E2F targets and G2M checkpoints as well as candidate drivers in cis contributing to enriched gene sets (*CENPF*, *ATF3)* (Fig. S8)*.* Together these observations imply an association between 1q+ and cell cycle dysregulation whether by enabling progression and entry into G2/M or inhibiting re-entry into G1. Driver genes involved in such cell cycle changes may include genes upregulated in G2 within 1q+ regions including *NUF2, NEK2*. 1q+ tumor subclones are also enriched in heat shock cell states and depleted in hypoxic cell states. Our integrated analysis of 1q+ showed a lower proportion of cells in hypoxia states within clones with 1q+ versus antecedent clones whereas the proportion of cells in a heat shock state was similar between clones (Fig. 4). The most abundant G1 state in clones with 1q+ showed elevated mitochondrial gene expression while cells in the hypoxia state showed notably decreased expression of mitochondrial genes.

We observed depletion of 16q- clones in post-mitotic states and enrichment in hypoxic states suggesting that 16q- may encourage cell cycle progression of post-mitotic cells into G1 in association with increased hypoxia-related gene expression. The diminished numbers of post-mitotic cells could also reflect a 16q- encouraged cell cycle arrest in G2 or M, although this is unlikely to be selected during tumorigenesis (Fig. 5). Differential expression comparisons between tumor subclones with and without 16q- revealed downregulation of metallothionein genes including *MT1X* and *MT2A*) not previously identified as candidate SCNA drivers (Fig. S10). Metallothionein synthesis can increase several-fold in response to oxidative stress, serving to safeguard cells against cytotoxicity, radiation, and DNA damage but this response is also associated with cell cycle exit (Ruttkay-Nedecky et al. 2013). Interestingly, because of their high cysteine content, MTs can readily sequester reactive oxygen species and protect cells against oxidative stress. It is suggested that these proteins, as a secondary antioxidant, cooperate with glutathione in maintaining the cellular redox state (Jamrozik et al. 2023).

In samples with 2p+ we found the formation of alternative cell states in G1 that clustered similarly in cell cycle space but showed divergent clone proportions. A G1 state that was depleted in 2p+ cells (here termed G1 A2) was marked by genes involved in cone photoreceptor differentiation including *PDC* showed downregulation of MYC targets and decreased expression of retinoblastoma molecular subtype 2 marker gene *TFF1* relative to the similarly positioned and 2p+ enriched G1A1(Fig. S14). This cluster was significantly depleted in the 2p+ clone in both tumors suggesting that avoidance of this state is a major 2p+ effect.

We also observed patterns of corresponding clusters enriched or depleted for 6p+ SCNAs with enrichment of progression-related expression. However, we did not observe common trends in gene set enrichment in the affected clusters according to 6p+ SCNA status. This may be due to variation in the 6p+ SCNA boundaries among our analyzed samples (Fig. S19, Table S5b) as well as differences in the co-occurring SCNAs or other genomic changes. Differential expression between clones across entire samples supports designation of previously proposed 6p+ candidate drivers including *DEK, E2F3, and SOX4* (Fig. S22), yet this pan-tumor upregulation could also indirectly result from other proliferation-inducing gene expression changes. Given findings of prognostic significance for 6p+ in ocular salvage (Berry et al. 2018) more scRNA-seq studies are urgently needed to power studies of 6p+ subclonal heterogeneity.

### SCNAs impact the formation of transcriptomic states

We identified cell states corresponding to cell cycle phases as well as stress-related states including hypoxia and hsp and observed variation in clone distribution according to cell state. To investigate continuous variation in expression programs characterizing cell states we used consensus non-negative matrix factorization (cNMF) with the mosaicMPI software package (Verhey et al. 2023). To identify factors that may explain effects of individual RB SCNAs we merged expression matrices from samples with specific copy number changes for each of 1q+, 2p+, 6p+, and 16q- (Fig. S23). This analysis showed close overlap between derived factors, cell cycle phase, and cell stress signals. Taken together, these clone distribution patterns suggest a scenario in which non-specific SCNAs induce stress states (mitochondrial or HSP expression) which are accommodated by RB SCNAs that allow cells to proceed through cell cycle checkpoints.. The evolutionary advantage conferred by SCNAs for tumor cells may be counterbalanced by the fitness penalty, such as proteasomal stress, that is expected to be associated with simultaneous dysregulation of many genes located on an affected region (Santaguida & Amon 2015). Given this fitness penalty, recurrent SCNAs likely include one or more driver genes that counterbalance strong negative selective pressures to enable cancer initiation, emergence of a malignant clone, or progression to more aggressive phenotypes.

We sought to identify whether SCNA clones might induce entirely new transcriptomic states. To identify clones that define an entirely new transcriptomic state we looked for those clones that made up a majority of a cluster (>70%) (Table S5). We found that the contribution of any single clone to a transcriptomic cluster rarely varied from the contribution of that clone to the tumor and did not exceed 70% (Table S6). Few transcriptomic clusters constituted entirely distinct SCNA clones. Instead, SCNA clones likely impact the distribution of cells across transcriptomic states in a probabilistic way.

Retinoblastoma molecular subtypes have been characterized by expression of later-stage cone differentiation markers (*ARR3*, *GUCA1C*) in one subtype and expression neuronal/ganglion cell markers (*SOX11, DCX*) (Liu et al. 2021). scRNA-seq analysis of a single retinoblastoma tumor suggests the existence of subclonal heterogeneity of molecular subtypes characterized by photoreceptor (subtype 1) or neuronal/ganglion expression (subtype 2) with distinct SCNAs (Liu et al. 2021). In our analysis, cell cycle identity strongly correlates with SCNA subclone identity in many tumors. Likewise, we observe a relationship between a subtype 2 gene signature score, aneuploidy, and cell cycle progression.

We next sought to identify SCNAs whose prevalence varied across transcriptomic clusters, as evidence that SCNAs induce alterations in cell state. We undertook a standard transcriptomic analysis of each tumor sample, starting with Louvain clustering at a range of resolutions. At low resolutions, clusters were identical to major cell cycle phases G1, S, G2 along with outlier clusters defined most of all by elevated expression of markers of heat shock and hypoxia. At higher resolutions, common intermediate cell cycle phases emerged recurrently in separate tumor samples. These more specific clusters often revealed finer points of cell cycle variation and were also enriched for specific SCNA clones. We observed that the most highly evolved tumor subclones appeared at highest frequency in proliferating transcriptomic states and in S/G2M cell cycle phases.

### SCNA subclones may exhibit altered expression within transcriptomic states

To identify transcriptomic changes associated with the acquisition of SCNAs, we differentiate between gene expression changes for coincident genes (in cis) and broader transcriptomic changes outside of the affected segment (in trans). Gene dosage effects in cis are likely to closely correlate with the degree of copy number change, though candidate driver genes may exhibit disproportionately large fold changes. Conversely, changes in trans may be indirectly influenced to a greater or lesser extent. Gene expression changes between two clones that differ solely by a specific SCNA can be attributed directly to SCNAs in the case of cis effects and indirectly in the case of trans effects. Numbers of differentially expressed genes *in cis* corresponded to the scale of the SCNA, with larger 1q+ and 6p+ regions showing a higher number of alterations. Interestingly differentially expressed genes *in trans* arose to the greatest extent in comparison between clones with and without 16q-. This finding supports prior designations of 1q+ and 16q- as distinctive features of aggressive retinoblastomas.

### Effects regulating G1 states

Across several tumor samples we identified distinct G1 transcriptomic states which frequently varied by clone composition. After integration of samples with 1q+ we observed enrichment of 1q+ clones in a dedifferentiated transcriptomic state G1 MT with elevated expression of mitochondrial genes and decreased expression of common photoreceptor genes including *GNB3 and DCT* (Fig. 4).

Direct comparison of G1 states showed enrichment of hallmark gene sets including *E2F targets* and *G2M checkpoints* and oncogenic signature gene sets including *RB P107 DN.V1 UP* in dedifferentiated states (Fig. S9) with elevated expression of *CENPF* and *CKS1B* in cis as well as *HCN3*. Notably, *HCN3* cooperates with *MYC* to promote mouse hepatocellular carcinoma and is frequently upregulated in human hepatocellular carcinoma where it is also co-increased with *MYCN* (Zhang et al. 2023). S state cells had enrichment of gene sets including APOPTOSIS and HYPOXIA, with elevated expression of differentially expressed genes including *ATF3, LMNA* (Fig. S8, Fig. S9).. It is possible that alterations in G1 gene expression may impact progression through later cell cycle stages. As an example, where different G1 states were characterized, Spencer et al. showed that as cells exit mitosis, they bifurcate into two populations with different levels of CDK2 activity (Spencer et al. 2013). Some cells rapidly increase their CDK2 activity and enter the next cell cycle immediately, while others lack CDK2 activity and go into a transient quiescent state according to CDK inhibitor p21 regulated during a restriction window at the end of the previous cell cycle when mitogens are present. If cells undergo stress including mitotic stress prior to G2 they will have prolonged G1 and otherwise will have short G1 (Spencer et al. 2013). Although specific variation in CDK activity is not clear in this work due to the detection limitations of scRNA-seq, elevated expression of cell cycle factors during G1 state is correlated with acquisition of 1q+ with upregulation of MYC targets gene sets and downregulation of differentiation-related expression.

### Candidate driver genes identified

This work attempts, for the first time, to identify retinoblastoma SCNA candidate driver genes by direct comparison of tumor subclones with and without specific SCNAs. Previous studies proposed candidate driver genes based on identification of minimal regions of gain and loss in array-based expression analysis and selection of target genes based on prior designation as cancer-related factors whether specifically in retinoblastoma or in cancer more generally (Herzog et al. 2001; Chen et al. 2001; Lillington et al. 2003; Zielinski et al. 2005; Gratias et al. 2007; Sampieri et al. 2009; Mol et al. 2014). A follow-up study validated a small number of candidate drivers through measurement of gene-dosage effects including *CRB1, NEK7* (1q), *MYCN* (2p), *SOX4, DEK* (6p) and multiple genes for 16q (Kooi et al. 2016a). We compared the set of differentially expressed genes between subclones with putative retinoblastoma SCNA driver genes, focusing on affected genes in cis (i.e., in the gained or lost segment) (Fig. 6). We similarly identified *DEK and SOX4 as well as CAP2, CCND3, E2F3, HLA-A, KIF13A, NUP153, PIM1, and TDP2* in 6p regionsas well as *KIF14, MCL1, MDM4, ZNF281* in 1q. It is important to note that the expression of proposed candidate driver genes may be correlated with cell cycle variation particularly *KIF14*, a Kinesin family member critical for regulation of p27 degradation and cell cycle regulation (Xu et al. 2014a). Finally, we identified several promising candidate driver genes not previously proposed. In 1q we found differential expression of *CENPF*, *ATF3* contributing toupregulation of hallmark gene sets including E2F targets and G2M checkpoint. In 16q we found downregulation of metallothionein genes (*MT1E, MT1X, MT1F, MT3, MT2A*) not previously identified as candidate SCNA drivers. In 2p we observed no significant differential expression of the frequently focally amplified *MYCN*. This omission may be due to limitations of detection within our scRNA-seq data where *MYCN* displays low (<0.5/cell) mean expression in all samples.

# Conclusions

In this work we sought to identify transcriptomic effects of somatic copy number alterations (SCNAs) in untreated retinoblastoma tumor scRNA-seq data. We identified and characterized tumor subclones with recurrent ‘RB SCNAs’ (1q+, 2p+, 6p+, and 16q-). We present a computational framework for analyzing tumor cell states and SCNA-associated transcriptional changes using Numbat to infer copy number alterations from droplet-based single-cell RNA sequencing data. Our analysis revealed that distinct cell clusters predominantly occupied specific positions within the cell cycle, expressing photoreceptor and neuronal markers during G1 phase, with phase-specific signatures during S and G2, along with stress-response markers including hypoxia and heat shock. We investigated potential relationships between SCNA occurrence and cell states, examining both the presence of SCNAs across different cell clusters and their relative enrichment or depletion within specific cell states compared to overall tumor frequencies. These findings demonstrated how Numbat-based analyses can illuminate distinct behavioral patterns of subclones within individual tumors. We detailed our comprehensive analysis of SCNA-associated transcriptional effects across all publicly available retinoblastoma single-cell RNA sequencing datasets. We examined 30 treatment-naïve retinoblastoma samples from five published studies (Collin et al. 2021; Yang et al. 2021; Wu et al. 2022; Field et al. 2022; Liu et al. 2024) using Numbat to infer copy number alterations. Our investigation characterized transcriptome profiles of retinoblastoma subclones with distinct SCNAs through analysis of tumor cell states, subclone-specific differential expression patterns, and identification of potential driver genes acting both in cis and trans.

Analysis of integrated scRNA-seq data from tumors with isolated chromosome 1q gain revealed a previously unidentified G2/M cell state predominantly composed of 1q+ subclones, along with enrichment of these 1q+ populations in post-mitotic and S/G2 cell cycle states (Fig. 4). Comparing gene expression between 1q+ and non-1q+ subclones within specific cell states highlighted upregulation of key gene sets including E2F targets and G2M checkpoints, as well as potential cis-acting driver genes such as CENPF and ATF3 that contribute to these enriched pathways (Fig. S8). These findings suggest that 1q gain influences cell cycle regulation, potentially by promoting G2/M progression or interfering with G1 re-entry. Analysis of 16q- subclones revealed their underrepresentation in post-mitotic states and overrepresentation in hypoxic states, suggesting that loss of 16q may promote cell cycle re-entry of post-mitotic cells into G1, coupled with enhanced hypoxia-related gene expression. Comparison of 16q- and non-16q- subclones identified novel candidate drivers: several downregulated metallothionein genes (MT1E, MT1X, MT1F, MT3, MT2A) not previously associated with SCNA effects (Fig. S10). In 2p+ samples, we identified distinct G1 cell states occupying similar cell cycle space but showing different clonal compositions. One G1 state, characterized by cone photoreceptor differentiation markers (GUK1, RXRG, SYT1, PDC) and reduced MYC target expression, showed consistent depletion of 2p+ cells across tumors, suggesting avoidance of this state is a key consequence of 2p gain. Notably, we observed no significant expression changes in the frequently amplified MYCN gene, possibly due to its low baseline expression (<0.5/cell) in our scRNA-seq data. While we identified clusters showing differential enrichment of 6p+ subclones and associated progression-related expression patterns, gene set enrichment analysis revealed no consistent trends across samples with 6p+. This variability may reflect differences in 6p+ SCNA boundaries among samples (Fig. S19, Table S5b). Given the prognostic importance of 6p+ in ocular salvage (Berry et al. 2018), additional scRNA-seq studies are needed to fully characterize 6p+ subclonal heterogeneity.

Together this work provides insight into the nature and extent of RB SCNA transcriptomic associations. The analysis of genetic subpopulations from single-cell RNA sequencing data provides an opportunity to simultaneously examine genetic and gene expression diversity as tumors evolve over time. Specifically, newly acquired changes in gene copy numbers can serve as innate genetic markers. When combined with distinctive patterns of gene expression, these markers enable researchers to monitor the dynamics of clonal groups throughout the course of tumor development. The current study is the first attempt in retinoblastoma to define transcriptomic effects of specific copy number changes by scRNA-seq analysis of tumor subclones. SCNAs in retinoblastoma are correlated with tumor progression with implications for therapeutic response, risk of metastatic disease, and survival. Identification of SCNA transcriptomic effects may have dramatic consequences for patients with other types of cancers. Our work demonstrates that precise characterization of SCNA effects by single cell analysis is possible given sufficient samples with subclonal heterogeneity. Further study is warranted to further elucidate the impact of RB SCNAs in retinoblastoma.

Limitations:

Do not account for SNVs. Tumors examined for epression of known variants (e.g., BCOR) but would not find RNAs degraded by nonsense mediated decay. Still, the most interpretable SCNAs were in tumors with limited SCNA evolution, likely also having least SNVs given correlation of exome SNVs and SCNAs in past datasets (may have to generate correlations based on accrued Stachelek datasets or see if others already did this).

Need more samples for more robust definititions, esp with SCNAs uniquely acquired over different pre-existing SCNA and SNV backgrounds. which will likely acrue over time.

Does not directly address different microenvironment (though if TEM has effect on cell state the alternative states should be seen and will be unique to samples in which the microenvironment effect is seen).

Still, despite limitations, results were robust enouch to enable discovery of differnetial cell cycle distributions, new states, and candidate driver genes.

Exclude previous candidate drivers (at least as drivers in these particular samples).

Depends on what tumors get sequenced, by chance.

# Materials and Methods

### Analysis of RB SCNAs in TCGA samples

The proportions of each cancer type with an RB SCNAs (1q+, 2p+, 6p+, and 16q-) were based on Table S2 of of Tayler et al (2018)(ref) and were defined as the sum of the relevant values (+1 for 1q+, 2p+, 6p+ and -1 for 16q-), divided by the number of samples of each type. GISTIC plots for 36 TCGA cancer types were produced by … .

### scRNA-seq data acquisition and quality control

scRNA-seq count matrices were retrieved for 31 samples from five studies (accession IDs GSE166175, PRJNA737188, GSE168434, GSE196420, GSE249995). fastq files were collected using SRA toolkit fasterq-dump. Cells were excluded from downstream Numbat and clustering analysis if percent mitochondrial reads was >5%, read counts <1,000, or genes detected < 1,000.

### Removal of non-tumor cells and clustering

Cells were clustered using the clustifyr R package (Swamy et al. 2021) and non-tumor cell clusters were identified and removed based on their diploid chromosome numbers assigned using Numbat and their cell types assigned using marker genes derived from the single cell Eye in a Disk database (ref).. After filtering cells were re-clustered with clustifyr, which gave different cluster assignments and marker genes differed from the unfiltered datasets (Fig. S4)..

### Numbat snakemake workflow

A complete snakemake workflow (Mölder et al. 2021) for SCNA inference and transcriptomic analysis is available at <https://github.com/whtns/external_rb_scrnaseq>. Pipeline software included Cell Ranger for alignment and read count quantitation, as well as Cellsnp-lite, eagle2 and samtools with respect to two different population genome reference panels: TOPMed and 1000 Genomes (1000G) as recommended by Numbat documentation for pileup and phasing (Auton et al. 2015; Loh et al. 2016; Taliun et al. 2021). We retrieved SNP VCF and phasing reference panels according to Numbat documentation. Cell-level phased allele counts were used for downstream Numbat inference along with per-sample count matrices from Cell Ranger processing.

To model expression profiles of contaminating normal retinal cell types in retinoblastoma scRNA-seq samples we used a fetal retinal scRNAseq dataset (Sridhar et al. 2020) as a reference control during SCNA inference. We evaluated the fitness of our chosen normal reference dataset by comparing expression noise score (mean squared error) using either fetal retinal cell types or using only cone cells. With noise minimized based on a reference containing all annotated retinal cell types. We ran Numbat for each sample with two iterations using default parameters except as noted below.

### Tuning Numbat parameters for successful SCNA inference

We iteratively tuned Numbat analysis using a reproducible snakemake workflow, with sample-specific parameter settings that were adjusted to accommodate each sample’s features and to recover likely SCNAs that were suggested by smoothed gene expression plots:

Initial clone number estimate (init\_k). Numbat uses hierarchical clustering of smoothed gene expression values to approximate an initial phylogeny given an initial clone number estimate (init\_k). To identify appropriate initial clone estimates for each sample, we consulted initial SCNA profiles inferred solely from smoothed gene expression data using a minimal Numbat model similar to InferCNV and obtained by … . When clones that were evident from visual inspection of smoothed gene expressio plots were not detected using the default init\_k: 3, we increased init\_k to aid detection of rare clones (*e.g.*, increasing from init\_k: 3 to init\_k: 4).

The minimum log-likelihood ratio (min\_LLR) is used to retain only high-confidence SCNAs in phylogeny reconstruction (default min\_LLR: 5). Here, we adjusted min\_LLR to 2 to recover lower-confidence CNVs when such SCNAs were visible in smoothed gene expression plots.

The maximum entropy parameter (max\_entropy) is used to filter SCNAs before phylogeny construction based on the uncertainty of an SCNA across individual cells. The entropy value ranges from 0 to 1, with 1 being the least stringent filter. We adjusted the original default max\_entropy: (0.5) to 0.8 to recover SCNAs that were visible in smoothed gene expression plots.

The formal name [*referred to as*] (Tau) specifies the stringency to simplify the mutational history (default tau: 0.3), which was adjusted upward to a maximum of ?? in order to recover SCNAs evident in smoothed gene expression plots but filtered out with default settings.

The transition probability (t) is used in the hidden Markov model to limit false-positive SCNAs, where a lower t is more appropriate for tumors with more complex copy number landscapes and a higher t is more effective for controlling false-positive rates of SCNA calls. The default t: 1e-5 was adjusted down to 1e-2 , as needed and in a sample-specific manner to recover SCNAs evident in smoothed gene expression plots.

**Inference of Phylogenies**

**Transcriptomic analysis**

We used Seurat (V4) (Butler et al. 2018; Stuart et al. 2019) to perform Louvain clustering and UMAP visualization of each sample. We performed preprocessing and Louvain clustering of retinoblastoma tumor scRNAseq datasets to identify gene expression-based cell clusters. We annotated single cell transcriptomes for each tumor sample through a common snakemake workflow including Cell Ranger for alignment and read count quantitation, Cellsnp-lite, eagle2 and samtools for pileup and phasing (Li et al. 2009; Huang & Huang 2021). Output from alignment and phasing were used for Numbat inference while quantitated UMIs were used for Seurat analysis. We scored cells using established cell-cycle phase marker genes (Tirosh et al. 2016) aggregated with Seurat *AddModuleScore* into S and G2/M gene signature scores. We found cluster marker genes by Wilcoxon rank-sum test and annotated common cell types likely arising due to normal cell contamination.

### Matrix factorization using mosaicMPI

To identify transcriptional states recurrent in retinoblastoma we employed consensus non-negative matrix factorization (cNMF) using the mosaicMPI software package (Verhey et al. 2023). For every individual sample we performed factorization at a range of k values from 6 to 24. To identify factors that may explain effects of individual RB SCNAs we merged expression matrices from samples with specific copy number changes for each of 1q+, 2p+, 6p+, and 16q-. Factorization of tumor samples is provided (Fig. S23).

### Regression of cell cycle effects

We compared transcriptomic clusters before and after regressing cell cycle effects. Regression was performed using Seurat CellCycleScoring and ScaleData function with phase markers according to the analysis from (Tirosh et al. 2016). We attempted to eliminate cell cycle phase effects before conducting clustering to enhance the visibility of recurrent cell states by removing what might be considered unimportant variation caused by the cell cycle. Correlation between cell cycle activity and cell states has been previously identified (Roccio et al. 2013; Soufi & Dalton 2016; Richard et al. 2018). We observed a complex interplay between cell cycle, clonal identity and tumor cell state, making cell cycle expression removal potentially problematic for accurate analysis.

### Comparison of clonal gene expression

We ran differential expression with MAST. We considered differentially expressed genes in cis with thresholds of adjusted p-value <= 0.05, log2 fold change >= 0.5 in the same direction as the SCNA of interest (Finak et al. 2015). We ran enrichment analyses of differentially expressed genes in cis and trans with clusterProfiler (Yu et al. 2012) using msigdb hallmark (H) and oncogenic gene sets (C6) and considering enriched terms with threshold adjusted p value <= 0.1.

# Acknowledgments

James Hicks

# References

Paste main manuscript references here in order of citation.

Automatic citation updates are disabled. To see the bibliography, click Refresh in the Zotero tab.

# Figures and Tables

***![](data:image/png;base64...)***

***Figure S1: RB SCNA frequency in TCGA by cancer type. (Taylor et al. 2018).***

**![](data:image/png;base64...)**

***Figure 2: Analysis of Numbat inference, cell state identity and subclone state distribution in a single tumor.***

*(A) SCNAs inferred with Numbat. (B) Clone inheritance trees. (C) UMAP plots colored by cluster, inferred cell cycle phase (Seurat), and SCNA subclone. (D) Cell cycle score plots colored by cluster, phase, and SCNA. (E) Faceted cell cycle score plots colored by SCNA. (F) Clone distribution plot colored by SCNA. Asterisks indicate significant differences from the tumor as a whole. \*, p<0.05 \*\*, p<0.01 \*\*\*, p<0.001, two-tailed binomial test.*

![](data:image/png;base64...)

***Figure 3: Analysis of cluster marker genes and cell state in a single tumor.***

*(A) Heatmap with tumor clusters and top 5 marker genes by log2FC. Annotated cluster identities on right and top of heatmap. (B) Tumor phylogeny by subclone inferred with Numbat. (C) Cells scored by gene signatures for S and G2M scores colored by cluster in cell-cycle space. (D) Cells plotted in cell cycle space faceted by cluster identity and colored by SCNA subclone. Cluster centroid labeled with colored diamonds.*

![](data:image/png;base64...)**Figure 4: 1q+ subclones in a multi-sample integration show redistribution to states G2/M and formation of novel G2/M state.**

*(A) Heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) Cell cycle score plots colored by cluster. (C) Faceted cell cycle score plots colored by SCNA subclone. (D) Clone distribution plot colored by SCNA subclone. \*, p<0.05; \*\*, p<0.01; \*\*\*, p<0.001; binomial test (E) Cell cycle distributions of subclones with 1q+ compared to antecedent clones in the integrated dataset (left) and individual tumors. Centroids indicate relative cluster sizes and positions. Tumor-specific phylogenies and karyograms are above and proportions of cells of each subclone in each cluster are below cell cycle plots, centroid and cluster colors as in panel A.*

![](data:image/png;base64...)

***Figure 5: 16q- subclones in a multi-sample integration show depletion in post-mitotic states and increased hypoxia.***

*(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster. (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone. (E) Proportion of cluster present in in subclones with 16q- compared to antecedent clones.*

![](data:image/png;base64...)

***Figure 6: RB SCNA candidate driver genes from prior studies.***

*Reproduced with modifications from Kooi et al. Red indicates clonal or subclonal gain. Blue indicates clonal or subclonal loss. Gray indicates proposed candidate drivers from prior studies. Yellow indicates significantly differentially expressed genes between clones with or without the relevant SCNA with concordant gene dosage changes. Average normalized expression across scRNA-seq samples displayed in green on bottom row.*

# List of Tables

**Table S1: RB SCNA frequencies in TCGA by cancer types.** Data extracted from Table S2 of Taylor et al. 2018.

**Table S2. Characteristics of retinoblastoma scRNA-seq datasets.**

**Table X2: Percent of clone in sample clusters.** (A) 1q+ (B) 16q- (C) 2p+. Significance of clone enrichment calculated by binomial test with respect to individual tumor and integrated analyses.

**Table S3: Rod cell proportions in tumor samples**

**Table S4: Tumor Quality Control.** Sequence read archive (SRA) ID, presence of subclonal heterogeneity, Number of UMIs, Number of genes detected, percent mitochondrial reads

**Table S5: Tumor metadata compiled from sequence read archive.** Includes sample ids, technical sequencing details, experimental metadata, and designation of all samples as untreated retinoblastoma.

**Table S6: SCNA segmentation from all Numbat analyses.** Chromosome boundaries and SCNA state (gain=amp, balanced gain=bamp, loss=del, loss of heterozygosity=loh)

**Table S7: Differentially expressed genes in in affected samples.** (A) within affected SCNA segments *in cis* (B) outside of affected segments *in trans*

**Table S8: Tally of removed clusters with marker genes.** Indicates removed diploid cell cluster sizes, marker genes, and cell type assignments for each tumor

**Table S9: Tally of candidate driver genes by differential expression of affected clone to antecedent clone.**

# Supplemental Methods

# Supplemental Data

**Supplementary Figure S1: Frequency of RB SCNAs in the TCGA dataset of 33 human cancers** . Data extracted from TCGA dataset (provide link).

**Supplementary Figure S2: Enrichment of RB SCNA segments in the TCGA dataset of 33 human cancers.** GISTIC plots of 33 human cancers with grey highlight of 1q, 2p, 6p, and 16q arms with significantly enriched copy number gain (red) or loss (blue) regions (q-values <??).

**Supplementary Figure S1: Quality control of scRNA-seq samples.**

Histograms display the number of UMIs and genes detected per cell and the percent of reads assigned to mitochondrial genes for all cells, color coded according to study. We retained 14 samples (including one technical replicate) from three studies and excluded 10 samples based on poor QC and seven samples (including two technical replicates) due to a lack of informative SCNAs. Technical replicates noted in brackets.

**Supplementary Figure S2: Expression-based SCNA profiles and hierarchical clustering of 11 retinoblastomas inferred by Numbat.**

Each row is an individual cell. Increased expression (relative to inferred diploid? cell) in red, decreased expression in blue. Gray horizontal lines indicate initial clone number parameter during Numbat execution.

**Supplementary Figure S3: Numbat SCNA inference heatmaps.**

SCNA subclones arranged vertically. Gains in red, losses in blue. Copy-neutral loss of heterozygosity (LOH) in green. Numbat inferred tumor clone phylogenies. Tumor samples arranged by presence of ‘RB SCNA’. Size of node indicates proportional cell number in subclone. Clone order inferred with Numbat. Karyograms with Numbat inferred RB SCNAs. Numbat inferred SCNA probability.

**Supplementary Figure S4: Removed non-tumor cell clusters identified by marker genes concordant with single cell Eye in a Disk database and the clustifyr R package.**

(A) Plots of cluster marker genes Rod cluster 5 (arrows) excluded in filtered dataset. (B) Paired UMAPS colored by cluster and clone assignment, with removed clusters identified (arrow).

**Supplementary Figure S5: Sample-specific analyses of tumors with 1q+ subclones after integration.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone.

**Supplementary Figure S6: Sample-specific analyses of tumors with 1q+ subclones without integration.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) faceted cell cycle score plots colored by SCNA subclone (C) UMAP plots colored by SCNA subclone (D) UMAP plots colored by cluster identity (E) clone distribution plot colored by SCNA subclone. (F) subclone phylogeny (G) cell cycle score plots colored by cluster.

**Supplementary Figure S7: Alternative resolutions for integrated 1q+ analysis.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone. The correspondence between cluster identity and cell cycle state varied with clustering resolution. At lower resolution, post-mitotic and G1 states occupied a common cluster while at higher resolution, many additional G1 states were distinguished. We observed no significant variation in distribution of 1q+ clones at either alternative resolution apart from enrichment of S/G2 and G2/M phases.

**Supplementary Figure S8: Differential expression of 1q+ clusters.**

(A) Candidate drivers found by differential expression between clones with 1q+ within clusters. (B) Enrichment analysis (over-representation) after differential expression comparisons between integrated 1q clusters of interest: G1 v. G1; S v. S\*; S/G2 v. G2/M

**Supplementary Figure S9: Enriched terms after differential expression comparison of G1 cell states in integrated 1q+:** G1 MT v. G1 DIFF and G1 MT v. G1 SUB2

(A) msigdb hallmark terms. (B) msigdb oncogenic signatures.

**Supplementary Figure S10: 16q- candidate drivers found by differential expression between clones within clusters.**

**Supplementary Figure S11: Alternative resolutions for integrated 16q- analysis.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone.

**Supplementary Figure S12: Sample-specific analyses of tumors with 16q- subclones after integration.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone. (E) Cluster similarity by distance of centroid in PCA.

**Supplementary Figure S13: Sample-specific analyses of tumors with 16q- subclones without integration.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) faceted cell cycle score plots colored by SCNA subclone (C) UMAP plots colored by SCNA subclone (D) UMAP plots colored by cluster identity (E) clone distribution plot colored by SCNA subclone. (F) subclone phylogeny (G) cell cycle score plots colored by cluster.

**Supplementary Figure S14: 2p+ subclones in a multi-sample integration show depletion in a G1 cell state.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster. (C) clone distribution plot colored by SCNA subclone (D) faceted cell cycle score plots colored by SCNA subclone (E) Proportion of cluster present in in subclones with 2p+ compared to antecedent clones.

**Supplementary Figure S15: Alternative resolutions for integrated 2p+ analysis.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone.

**Supplementary Figure S16: Sample-specific analyses of tumors with 2p+ subclones after integration.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) cell cycle score plots colored by cluster (C) faceted cell cycle score plots colored by SCNA subclone (D) clone distribution plot colored by SCNA subclone. (E) Cluster similarity by distance of centroid in PCA.

**Supplementary Figure S17: Differential expression and enriched hallmark terms between 2p+ enriched g1 cluster versus 2p+ depleted clusters.**

(A-B) g1 a1 v. g1 a2 (C-D) g1 b v. g1 a2 (E-F) g1 b v. g1 a1

**Supplementary Figure S18: Candidate driver genes in 2p+ in cell states.**

**Supplementary Figure S19: 6p+ SCNA boundaries in affected samples.**

**Supplementary Figure S20: Sample-specific analyses of tumors with 6p+ subclones without integration.**

**Supplementary Figure S21: Enriched terms in comparisons between 6p+ enriched and 6p+ depleted clusters.**

In tumor SRX14116944 (with gain of all of 6p superimposed over 1q+), 6p+ and non-6p+ clones occupied distinct G1 states with *CRABP2* (a factor that mediates RXRg signaling) and subtype 2 retinoblastoma marker gene *CD24* characterizing 6p+ enriched clusters and cone differentiation markers *PDE6H and PDC* characterizing 6p+ depleted clusters (Fig. S20). Direct comparison of clusters within the same G1 cell cycle phase with and without enrichment of 6p+ clones showed activation of HALLMARK MYC TARGETS and suppression of TNFA SIGNALING VIA NFKB in G1 states in the 6p+-enriched cluster (Fig. S21).

In tumor SRX10264524 (with concurrent gain of 6p sequences from ~5 MB to ~ 38 MB (Chr6:5-38 MB) and 11p loss superimposed over diploid), 6p+ clones also occupied distinct G1 states with enrichment in clusters characterized by expression of neurogenic marker *GAP43* and depletion in clusters characterized by stemness marker *CD24* and mature cone markers including *GUCA1A and ARR3* (Fig. S20) Direct comparison of clusters within the same G1 cell cycle phase

with and without enrichment of the 6p+ in clones showed activation of E2F targets and suppression of TGF BETA SIGNALING in G1 states (Fig. S21). The very existence of variation in 6p+ clone distribution across cell states suggests an impact of the 6p region from 5-38 MB on the transcriptome, although transcriptomic changes could also relate to the concurrent 11p loss or undetected co-occurring mutations.

**Supplementary Figure S22: Differentially expressed genes in cis after comparisons between 6p+ subclones and antecedent clones across all clusters in tumor.**

Genes differentially expressed between initial 6p+ clones and antecedent non-6p+ clones in two retinoblastoma tumors and all genes positioned within the 6p-gain segment Chr6:0-69 MB (SRX14116944) and the Chr6:5-38 MB (SRX10264524).

**Supplemental Figure S20: Sample-specific analyses of tumors with 6p+ subclones without integration.**

(A) heatmap of cluster marker genes defined by Wilcoxon rank-sum test. (B) faceted cell cycle score plots colored by SCNA subclone (C) UMAP plots colored by SCNA subclone (D) UMAP plots colored by cluster identity (E) clone distribution plot colored by SCNA subclone. (F) subclone phylogeny (G) cell cycle score plots colored by cluster.

![A screenshot of a computer screen

Description automatically generated](data:image/png;base64...)![A screenshot of a computer

Description automatically generated](data:image/png;base64...)![A screenshot of a computer screen

Description automatically generated](data:image/png;base64...)![A screenshot of a computer screen

Description automatically generated](data:image/png;base64...)

**Supplemental Figure S23: cNMF factorization and clustering of samples with RB SCNAs.**

**Displayed samples containing with subclonal heterogeneity for RB SCNAs.** 1q: SRX10264526, SRX11133594, SRX11133593, SRX11133592; 2p+: SRX10264523, SRX10264524, SRX10264525, SRX14116944; 6p+: SRX10264524, SRX10264525, SRX14116944, SRX22868105; 16q-: SRX11133594, SRX11133593, SRX11133592. Displayed SCNAs include subclones with relevant SCNA and antecedent subclones (see Figure 4.2). Column annotation rows display factorization values for 6 factors ![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a graph

Description automatically generated](data:image/png;base64...)![A chart of a number of colors

Description automatically generated with medium confidence](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A chart of a number of gene

Description automatically generated with medium confidence](data:image/png;base64...)![A screen shot of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a graph

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)![A close-up of a chart

Description automatically generated](data:image/png;base64...)

![A screenshot of a computer

Description automatically generated](data:image/png;base64...)

***Figure S24: Numbat inferred tumor clone phylogenies.***

*Tumor samples arranged by presence of ‘RB SCNA’. Size of node indicates proportional cell number in subclone. Clone order inferred with Numbat. Excluded tumors were not evaluated for the effects of indicated SCNAs for the reasons indicated, although they may be evaluated for effects of other SCNAs* (*e.g., SRR17960484 evaluated for 2p+ and 6p+ but not 1q+). “Few preceding” = cells in antecedent clones insufficient for differential expression. “With confounding SCNA” = co-inherited RB SCNAs. “No clone distribution” = differential expression and clone distribution comparisons not powered for analysis*

![](data:image/png;base64...)

***Figure S25: Subtype violins***

*violin plot with subtype 2 gene signature scores for G1 clusters to bolster the distinction*

![](data:image/png;base64...)

***Figure S26: CENPF violins***

*violin plot with CENPF scores for 1q+ samples*

# Extra

### Analysis of SCNA effects on transcriptomic cell state distribution

We evaluated whether Numbat can be used to detect SCNA-associated changes in transcriptomic cell states or in the distribution of tumor cells between such states. In a demonstration analysis, we examined tumor SRX11133594 (Fig. 5) in which Numbat identified three subclones: a diploid subclone with no SCNAs, acquired 16q-, and acquired 16q- and 1q+ (Fig. 5A,B). Then, using Louvain clustering, we visualized relationships between clusters in UMAP plots (Becht et al. 2019).

UMAP plots segregated cells according to their Louvain clusters but showed limited segregation according to subclone composition (Fig. 5C), suggesting that variation due to subclone identity may be secondary compared to other gene expression determinants, such as cell cycle phase. To assess the relationship of cell cycle phase to subclone composition, we colored the UMAP plot according to each cell’s G1, S, and G2/M assignment using the average expression levels of gene sets associated with cell cycle phases subtracted by the aggregated expression of control gene sets. (Tirosh et al. 2016). However, simply categorizing cells into phase based on simplified assigned cell cycle phase (G1, S, G2) resulted in assignment to cell cycle phases that elided distinct variations according to cluster (Fig. 5D). This may be because assigning cells to only three cell cycle states is somewhat arbitrary; although cell cycle position is an important determinant it is not as simple as G1, S and G2/M (Guo & Chen 2024). In contrast, plotting the different Louvain clusters in cell cycle space enabled comparisons according to cell cycle position (Fig. 5D). Because clusters overlap, we marked average cell cycle position with centroids and separately displayed clusters to provide an overview of each cluster’s cell cycle distribution (Fig. 5D, E).

We colored cells according to their SCNA subclone status to visually display the different clone distributions (Fig. 5E) and quantified the clone distributions in each state with stacked bar plots, where distributions are shown relative to the distributions of all tumor subclones of interest combined (Fig. 5F). This allowed us to address whether acquisition of RB SCNAs was associated with emergence of new transcriptomic states or with an altered distribution of tumor cells (Fig. 5F).For this tumor, it was evident that both occurred, as each cluster had a distinct subclone composition and a G2/M cluster was evident solely in the 16q-.1q+ subclone. For most clusters, the subclone distributions significantly differed from the tumor as a whole based on a binomial test.

In analysis of this and other individual tumors, we identified cell clusters centered in G0/G1, G1/S, S/G2, G2/M, and post-mitotic phases. However, some clusters, such as the two G1/S clusters in Fig. 5E, had similar cell cycle distributions and centroid positions, but different subclone distributions, suggesting that there were two versions of the G1/S state that were differentially populated according to the SCNA status.

![A screenshot of a computer screen

Description automatically generated](data:image/png;base64...)

***Figure 5: Analysis of Numbat inference, cell state identity and subclone state distribution in a single tumor.***

(A) SCNAs inferred with Numbat. (B) Clone inheritance trees. (C) UMAP plots colored by cluster, inferred cell cycle phase (Seurat), and SCNA subclone. (D) Cell cycle score plots colored by cluster, phase, and SCNA. (E) Faceted cell cycle score plots colored by SCNA. (F) Clone distribution plot colored by SCNA. Asterisks indicate significant differences from the tumor as a whole. \*, p \*\*, p<. \*\*\*, p< , two-tailed binomial test.

While the precise marker genes varied for different clusters in different tumors, we distinguished G0 and/or G1 states that differed according to expression of markers of cone photoreceptor differentiation such as *GNB3* (Ritchey et al. 2010) and *GUK1* (Kallman et al. 2020). We also distinguished two S phase clusters in most samples: a major S cluster often marked by high expression of S phase marker genes *HMGB1 (*high mobility group box 1*), RRM2* (ribonucleotide reductase), and various histone genes, and a smaller population distinguished by high expression of S and M-related chromatin regulators *UBE2T* and *ATAD2* and having a higher S phase score than the more abundant S phase cells, which we labeled S+. S/G2 clusters have additional G2-related gene expression including TOP2A, CENPF, and UBE2C. We also frequently identified a post-mitotic (pm) cell cluster characterized by expression of G2 signature with low S score and expression of residual G2 genes including *PTTG1* and *ARL6IP1*. We also identify transcriptomic states characterized by markers of heat shock (e.g., HSPA1A or HSPB1 (O’Flanagan et al. 2019)), and hypoxia (BNIP3, ENO1) with cell cycle positions spanning g1/s/g2 or centered on g1, respectively (Fig. 5). Thus, cell cycle variation was a major contribution to distinct cell states along with stress states distinguished by hypoxia and heat shock markers.

We therefore concluded that retinoblastoma cell states might be intrinsically linked to cell cycle variation, consistent with a previously identified correlation between cell cycle activity and cell states (Roccio et al. 2013; Soufi & Dalton 2016; Richard et al. 2018). Thus, we observed a complex interplay between cell cycle, clonal identity and tumor cell state, making cell cycle expression removal potentially problematic for accurate analysis. On the other hand, retaining cell cycle position information enabled interpretable comparisons in cell cycle distributions of the different SCNA clones.

### Tumor subclones with 6p+ are enriched in a G1 cell state with diverse expression profiles

We identified four potentially interpretable samples affected by 6p+ (Fig. S24). 6p+ was an initiating or co-initiating event in three tumors (SRX10264524, SRX10264525, SRX22868105) and was preceded by 1q+ in one case (SRX14116944). We examined boundaries of 6p+ changes in the four proposed samples and identified SCNA boundaries ranging from total 69 MB telomere-to-centromere bound (SRX14116944), 10 MB telomere-bound (SRX10264525), 33 MB unbound (SRX10264524), and 8 MB centromere-bound (SRX22868105) (Fig. S19, Table S5b). The varying 6p+ boundaries and variable preceding or co-occurring SCNAs complicated comparison between samples. In fact, integration of the 6p+ and antecedent clones from multiple samples did not reveal common patterns within 6p+ subclones across samples. Recent studies suggest that recurrent SCNAs are selected during tumorigenesis due to their altered expression of cancer driver genes rather than because of chromosome-biased rates of mitotic segregation errors (Shih et al. 2023). Therefore, the presence of non-overlapping 6p+ regions (e.g., in SRX10264525 and SRX22868105) may point to entirely different drivers and transcriptome consequences. We therefore examined tumors individually and found cell states corresponding to cell cycle phase, including G1, S, G2, G2/M, as well as stress states hypoxia and heat shock (Fig. S20). However, we observed variation in clone distribution according to cell state only in samples SRX14116944 and SRX10264524. The lack of variation in clone distribution in SRX22868105 and SRX10264525 may be due to the limited number of cells in the 6p+ and antecedent clones, which may have diminished the ability to segregate 6p+-enriched and 6p+-depleted versions of cell cycle states and the smaller SCNA boundaries <=10Mb which are uncharacteristic of typically arm-level 6p+ changes in retinoblastoma.

To identify candidate driver genes in cis we performed differential expression between clones with and without 6p+ for each cluster and for the different 6p-gain regions in each of the above samples. The comparison yielded significantly differentially expressed genes in cis,with several cell cycle related genes including *KIFC1* and dramatic variation in known 6p+ candidate driver *DEK*, but no significant upregulation of previously proposed candidates E2F3 and SOX4 *(*Fig. S22, Table S7b*).* Notably, some 6p genes outside the region gained in SRX10264524 were significantly upregulated in the SRX10264524 clone, suggesting these may be induced secondarily as a consequence of the gains in the Chr6:5-38 MB region.

These limited findings in two samples with overlapping 6p+ regions suggest a role of one or more 6p+ regions in generating new transcriptomic states characterized by increased stemness and neurogenic signaling, The decreased contribution of 6p+ clones to G1 states marked by cone photoreceptor terminal differentiation markers *PDE6H and PDC* (for SRX14116944) or *GUCA1A* (for SRX10264524) suggests that gain of the Chr6:5-38 MB region suppresses terminal differentiation to enable continued cell cycling. However, concrete findings are limited by lack of available samples for direct with-versus-without-6p+ clone comparisons due to frequent co-occurring, confounding SCNAs.

## Strategies to identify SCNA candidate driver genes *some of this section elaborates on points made in the first Intro paragraph, bawed on which I copy-pasted some of the references to that paragraph (though can you check if they’re cited correctly?). I also ‘moved up’ aspects more related to SCNA effects than to driver genes per se (which were for us a secondary focus). Can you check if other points would be helpful to include and otherwise delete w track changes?*

Candidate driver genes can be identified by minimal recurrence analysis, correlated sequencing data, and functional tests.

Specific aneuploidies clearly influence cancer initiation, progression, and resistance, but their precise roles are complex and context dependent.

Retinoblastoma is a pediatric retinal cancer that may be used to dissect the role of progression-related subclonal mutations, as most retinoblastomas initiate in response to the same genomic alteration, pass through similar pre-malignant stages, and carry low mutational burdens. Specifically, ~98% of retinoblastomas are thought to form in response to biallelic inactivation of *RB1* in maturing cone precursors and to proceed through a premalignant stage before retinoblastomas appear, whereas ~ 2% form in immature cone precursors in response to *MYCN* amplification (Rushlow et al. 2013; Xu et al. 2014b; Singh et al. 2018, 2022) *RB1*-mutant retinoblastomas may initially lack chromosomal abnormalities or single nucleotide variants (SNVs) beyond *RB1* (Kooi et al. 2016c), implying that progression-related genomic changes are superimposed over minimally altered genomic landscapes.


