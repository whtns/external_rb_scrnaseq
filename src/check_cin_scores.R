# seu <- readRDS("output/seurat/SRR14800534_unfiltered_seu.rds")

# seu <- readRDS("output/seurat/SRR13884247_unfiltered_seu.rds")

score_chrom_instability <- function(seu_path) {

  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)

  cin_scores <- list(
  "cin70" = c("TPX2", "PRC1", "FOXM1", "CDC2", "TGIF2", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "ch-TOG", "NUP205", "CDC20", "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "BRRN1", "MCM10", "LSM4", "MTB", "ASF1B", "ZWINT", "TOPK", "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A", "GPIandMGC13096", "SFRS2", "HDGF", "NXT1", "NEK2", "DHCR7", "STK6", "NDUFAB1", "KIAA0286", "KIF4A"
  ),
  "pos_tri70" = c("ANXA7", "ATG7", "BDNF", "CDKN2B", "CHFR", "CNDP2", "ENO3", "F3", "GIPC2", "GJB5", "HSPB7", "IMPACT", "P2RY14", "PCDH7", "PKD1", "PLCG2", "SNCA", "SNCG", "TMEM140", "TMEM40"),
  "neg_tri70" = c("AURKA", "BCL11B", "BIRC5", "BLMH", "BRD8", "BUB1B", "CCNA2", "CDC5L", "CDK1", "CENPE", "CENPN", "CTH", "DLGAP5", "HMGB2", "IDH2", "ISOC1", "KIAA0101", "KIF22", "LIG1", "LSM2", "MCM2", "MCM5", "MCM7", "MYBL2", "NASP", "NCAPD2", "NCAPH", "NMI", "NUDT21", "PCNA", "PIGO", "PLK1", "PLK4", "POLD2", "RACGAP1", "RAD51", "RFC2", "RFC3", "RFC5", "RPA2", "SMAD4", "SMC4", "SSRP1", "TAB2", "TCOF1", "TIMELESS", "TIPIN", "TOP2A", "UBE2C", "USP1"),
  "het70" = c("AHCYL1", "AKT3", "ANO10", "ANTXR1", "ATP6V0E1", "ATXN1", "B4GALT2", "BASP1", "BHLHE40", "BLVRA", "CALU", "CAP1", "CAST", "CAV1", "CLIC4", "CTSL1", "CYB5R3", "ELOVL1", "EMP3", "FKBP14", "FN1", "FST", "GNA12", "GOLT1B", "HECTD3", "HEG1", "HOMER3", "IGFBP3", "IL6ST", "ITCH", "LEPRE1", "LEPREL1", "LEPROT", "LGALS1", "LIMA1", "LPP", "MED8", "MMP2", "MUL1", "MYO10", "NAGK", "NR1D2", "NRIP3", "P4HA2", "PKIG", "PLOD2", "PMP22", "POFUT2", "POMGNT1", "PRKAR2A", "RAGE", "RHOC", "RRAGC", "SEC22B", "SERPINB8", "SPAG9", "SQSTM1", "TIMP2", "TMEM111", "TRIM16", "TRIO", "TUBB2A", "VEGFC", "VIM", "WASL", "YIPF5", "YKT6", "ZBTB38", "ZCCHC24", "ZMPSTE24"),
  "bucc_up" = c("HPDL", "L1CAM", "SMPD4", "ELF4", "IRAK1", "ISG20L2", "ZNF512B", "TMEM164", "RCE1", "AMPD2", "SCXB", "HSF1", "PLXNA3", "RBP4", "MRPL11", "NSUN5", "PNMA5", "FADS3", "TRAIP", "IGSF9", "TFAP4", "BOP1", "FAM131A", "EIF1AD", "FAM122C", "SLC10A3", "CCNK", "ELK1", "DGAT1", "DHCR7", "GNAZ", "SCAMP3", "ZNF7", "CCDC86", "SLC35A2", "CPSF4", "USP39", "FTSJD2", "DBF4B", "MEPCE", "ZBTB2", "SEMA4D", "LRRC14", "YKT6", "METTL7B", "SCAMP5", "IP6K1", "MTCH1", "SLC6A9", "SMCR7L", "SAMD4B", "CD3EAP", "PPP1R16A", "TMCO6", "U2AF2", "POM121", "TSR2", "UNKL", "IFRD2", "CYHR1", "MAGEA9B", "TBC1D25", "ZNF275", "TJAP1", "CLDN4", "SPNS1", "TMEM120A", "CPSF1", "FABP3", "FAM58A"),
  "bucc_down" = c("CYTIP", "HLA-DMB", "FYCO1", "NUBP1", "AGPS", "CADPS2", "C3orf52", "FGL2", "PDIA3P", "IL18", "APOBEC3G", "CA8", "HNRNPA1L2", "CD44", "TMEM2", "ELL2", "FCGBP", "TAPT1", "AHR", "HLCS", "AEN", "AGR3", "STEAP1", "TLR4", "HNRNPA1", "HLA-DRB5", "SLITRK6", "LPXN", "HLA-DRB1", "BACE2", "SYTL1", "KYNU", "RPS27L", "CD96", "HLA-DPA1", "TRA2A", "FBXO25", "HFE", "DRAM1", "DYRK4", "FBXL5", "GLUL", "HLA-DRA", "KCNMA1", "IL18R1", "SIPA1L2", "STEAP2", "TP53I3", "LIMA1", "BMX", "HLA-DMA", "C17orf97", "RPH3AL", "SLC37A1", "PHLDA3", "TK2", "NHS", "TMEM229B", "CYFIP1", "B3GALT5", "FAM107B", "ADH6", "CCND1", "LCK", "RABL2B", "SPA17", "GADD45A", "FHL2", "RPL36AL", "CCDC60", "MUC5B", "PHLDA1", "ZG16B", "SIDT1", "BTD", "TTC9", "GSTZ1", "ITGA6", "C10orf32", "FAM46C", "SPDEF", "ASRGL1", "TCN1", "NUDT2", "CDKN1A", "POLR3GL", "C4BPA", "LAMC2", "SLC6A14")
)

  seu <- AddModuleScore(seu, sheltzer_scores, name = "cin")

  names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("cin", seq(1, length(cin_scores))))] <- names(cin_scores)

  cin_score_fplot <- FeaturePlot(seu, names(cin_scores))

  plot_path <- glue("results/effect_of_regression/{sample_id}_cin_scores.pdf")

  ggsave(plot_path, height = 8, width = 8)

  return(plot_path)
}

score_chrom_instability("output/seurat/SRR13884242_unfiltered_seu.rds")

