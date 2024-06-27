source("packages.R")
source("functions.R")
library(targets)
library(org.Hs.eg.db)



seu <- readRDS("output/seurat/SRR13884246_branch_5_filtered_seu.rds")

seu$clusters

# g1 ------------------------------
g1_diffex <- Seurat::FindMarkers(seu, group.by = "clusters", ident.1 = "g1_0", ident.2 = "g1_2") %>% 
	clean_diffex()

down_csv(g1_diffex)

debug(enrichment_analysis)

# s diffex ------------------------------

table_and_plot_enrichment <- function(seu, groups) {
	
	my_diffex <- Seurat::FindMarkers(seu, group.by = "clusters", ident.1 = groups[[1]], ident.2 = groups[[2]])
	
	my_enrichment <- 
		my_diffex %>% 
		# clean_diffex() %>% 
		enrichment_analysis(gene_set = "hallmark") %>% 
		setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
	
	my_enrichment_plot <- 
		my_enrichment %>% 
		plot_enrichment() + 
		labs(title = glue("{names(groups)[[1]]}: {groups[[1]]} v. {names(groups)[[2]]}: {groups[[2]]}"))
	
	netplot <- cnetplot(my_enrichment, node_label="all", 
					 showCategory = 10, 
					 # cex.params = list(gene_node = 2), 
					 cex_category= 0.5,
					 cex_gene= 0.5,
					 cex_label_category= 0.5,
					 cex_label_gene= 0.5
	)
	
	return(
		list("diffex" = clean_diffex(my_diffex),
				 "enrichment_table" = my_enrichment, 
				 "dotplot" = my_enrichment_plot,
				 "netplot" = netplot
				 )
		)
}



effects_2p_gain_g1 <- table_and_plot_enrichment(seu, groups = list("2p+" = "g1_0",
																						 "diploid" = "g1_2"))
effects_2p_gain_g1$dotplot + 
	effects_2p_gain_g1$netplot + 
	plot_layout(widths = c(1,3))

ggsave("results/effects_2p_gain_g1_plot.pdf", height = 10, width = 16) %>% 
	browseURL()

effects_2p_gain_s <- table_and_plot_enrichment(seu, groups = list("2p+" = "s_g2_7",
																						 "diploid" = "s_6"))

effects_2p_gain_s$dotplot + 
	effects_2p_gain_s$netplot + 
	plot_layout(widths = c(1,3))

ggsave("results/effects_2p_gain_s_plot.pdf", height = 10, width = 16) %>% 
	browseURL()



