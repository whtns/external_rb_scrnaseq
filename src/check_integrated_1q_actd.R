library(biomaRt)

library(gprofiler2)

seu <- readRDS("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")

mouse_iegs <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3", "Egr4", "Nr4a1", "Nr4a2", "Nr4a3", "Arc", "Homer1", "Rheb", "Rgs2", "Plk2", "Ptgs2", "Bdnf", "Inhba", "Nptx2", "Plat", "Nrn1", "Myc", "Dusp1", "Dusp5", "Dusp6", "Pcdh8", "Cyr61", "Gadd45b", "Trib1", "Gem", "Btg2", "Ier2", "Npas4", "Rasd1", "Crem", "Mbnl2", "Arf4", "Gadd45g", "Arih1", "Nup98", "Ppp1r15a", "Fbxo33", "Per1", "Per2", "Maff", "Zfp36", "Srf", "Mcl1", "Ctgf", "Il6", "Atf3", "Rcan1", "Ncoa7", "Cxcl2", "Bhlhe40", "Slc2a3", "Nfkbia", "Ier3", "Sgk1", "Klf6", "Klf10", "Nfkbiz", "Flg", "Gbp2b", "Tnfaip3", "Cebpd", "Hbegf", "Ldlr", "Tsc22d1", "F3", "Ccl2", "Csrnp1", "Pmaip1", "Zfp36l2", "Plau", "Ccl5", "Saa3", "Ifnb1", "Tnf", "Irf1", "Cd83", "Map3k8", "Socs3", "Csf2", "Il1a", "Cxcl1", "Il12b", "Il1b", "Sod2", "Pim1", "Peli1", "Tlr2", "Ccl3", "Noct", "Bcl3", "Ifit2", "Icam1", "Ifit1", "Tnfsf9", "Ccrl2", "Cxcl10", "Gbp2", "Il10", "Clec4e", "Acod1", "Mmp13", "Cxcl11", "Il23a", "Arhgef3", "Serpine1", "Traf1", "Vcam1", "Ackr4", "Marcksl1", "Nfkbid", "Ikbke", "Ccl12", "Ifit3", "Cebpb", "Zfp36l1", "Txnip", "Nfib", "Hes1", "Pias1", "Klf2", "Cd69", "Dusp2", "Wee1", "Thbs1", "Sik1", "Gdf15", "Ier5", "Rgs1", "Id2", "Apold1")

human_iegs = gorth(mouse_iegs, source_organism = "mmusculus", target_organism = "hsapiens")$ortholog_name

DefaultAssay(seu) <- "gene"

seu <- Seurat::AddModuleScore(seu, list(human_iegs), name = "actd")

DefaultAssay(seu) <- "integrated"

plot_path <- "results/check_integrated_1q_actd.pdf"
pdf(plot_path)

VlnPlot(seu, features = "actd1", group.by = "integrated_snn_res.0.4")
VlnPlot(seu, features = "actd1", group.by = "integrated_snn_res.0.6")
VlnPlot(seu, features = "actd1", group.by = "integrated_snn_res.0.8")

dev.off()

browseURL(plot_path)
