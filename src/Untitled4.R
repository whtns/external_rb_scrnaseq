
library(targets)
source("packages.R")
source("functions.R")

tar_load(c("corresponding_clusters_diffex", "corresponding_seus"))

plot_corresponding_enrichment <- function(sample_diffex_list, sample_id, ...){
	# browser()
	h_gene_set = msigdbr(species = "human", category = "H") |> 
		dplyr::select(gs_name, entrez_gene)
		
	c6_gene_set = msigdbr(species = "human", category = "C6") |> 
		dplyr::select(gs_name, entrez_gene)

	c6_plots <- sample_diffex_list |> 
		map("all") |> 
		map(prep_for_enrichment, TERM2GENE = c6_gene_set) |> 
		map(plot_enrichment) |> 
		imap(~{.x + labs(title = sample_id, subtitle = .y, caption = "c6")})
	
	h_plots <- sample_diffex_list |> 
		map("all") |> 
		map(prep_for_enrichment, TERM2GENE = h_gene_set) |> 
		map(plot_enrichment) |> 
		imap(~{.x + labs(title = sample_id, subtitle = .y, caption = "h")})
	
	plot_path <- glue("results/{sample_id}_corresponding_clusters_diffex_enrichment.pdf", h = 10)
	pdf(plot_path)
	print(c6_plots)
	print(h_plots)
	dev.off()
	
	return(plot_path)
		
}

test0 <- 
corresponding_clusters_diffex %>%
	set_names(fs::path_file(unlist(corresponding_seus))) |> 
	imap(plot_corresponding_enrichment)


plot_corresponding_enrichment(corresponding_clusters_diffex$corresponding_clusters_diffex_15edb65d$`g1_1-g1_0-g1_3-g1_11 v. g1_2-g1_5`)


imap(corresponding_clusters_diffex$corresponding_clusters_diffex_15edb65d)



test0 <- map(corresponding_clusters_diffex, enrich_corresponding_clusters, TERM2GENE = c6_gene_set)



test0 <- map_depth(corresponding_clusters_diffex, 2, "all") |> 
	map(list_flatten) |> 
	map_depth(2, ~{plot_enrichment(prep_for_enrichment(.x, c6_gene_set))}) |>
	identity()

# myenrichment <- 
	corresponding_clusters_diffex$corresponding_clusters_diffex_15edb65d$`g1_1-g1_0-g1_3-g1_11 v. g1_2-g1_5`$all |> 
	prep_for_enrichment(TERM2GENE = c6_gene_set) |> 
	plot_enrichment()

# myenrichment <- 
	corresponding_clusters_diffex$corresponding_clusters_diffex_15edb65d$`s_g2_7-s_9 v. s_g2_8`$all |> 
	prep_for_enrichment(TERM2GENE = c6_gene_set) |> 
	plot_enrichment()

# myenrichment <- 
	corresponding_clusters_diffex$corresponding_clusters_diffex_5b95071f$`g1_3-g1_6 v. g1_5`$all |> 
	prep_for_enrichment(TERM2GENE = c6_gene_set) |> 
	plot_enrichment()

corresponding_clusters_diffex$corresponding_clusters_diffex_40baf62d$`g1_1 v. g1_6-g1_7`$all |> 
prep_for_enrichment(TERM2GENE = c6_gene_set) |> 
plot_enrichment()

corresponding_clusters_diffex$corresponding_clusters_diffex_a6e7c676$`g1_3 v. g1_1-g1_5`$all |> 
	prep_for_enrichment(TERM2GENE = c6_gene_set) |> 
	plot_enrichment()

corresponding_clusters_diffex$corresponding_clusters_diffex_a6e7c676$`s_6 v. s_7`$all |> 
	prep_for_enrichment(TERM2GENE = c6_gene_set) |> 
	plot_enrichment()
