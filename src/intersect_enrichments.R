
tar_load(c("corresponding_seus", "corresponding_clusters_enrichments"))

enrichments <- 
	map(corresponding_clusters_enrichments, "enrichment") |> 
	set_names(corresponding_seus) |> 
	identity()

enrichments_2p <-
	enrichments[2:4] |>  
	purrr::list_flatten()

enrichments_2p_c6 <- 
	enrichments_2p[str_detect(names(enrichments_2p), "_c6")] |> 
	purrr::list_flatten()

enrichments_2p_h <- 
	enrichments_2p[str_detect(names(enrichments_2p), "_h")] |> 
	purrr::list_flatten()

enrichments_2p_h <- enrichments_2p_h[str_detect(names(enrichments_2p_h), "g1")]

common_terms <- map(enrichments_2p_h, ~{.x@result[["ID"]]}) |> 
	purrr::compact()

common_terms_table <-
	common_terms |> 
	enframe("analysis", "terms") |> 
	dplyr::mutate(sample_id = str_extract(analysis, "SRR[0-9]*")) |> 
	dplyr::mutate(comparison = str_extract(analysis, "_g1.*v.*$")) |> 
	tidyr::unnest(terms) |> 
	dplyr::group_by(terms) |> 
	dplyr::mutate(recurrence = n_distinct(sample_id)) |> 
	dplyr::arrange(desc(recurrence), terms) |> 
	dplyr::filter(recurrence == 3) |> 
	dplyr::pull(terms) |> 
	unique()

filter_enrichment_by_term <- function(enrich_out, terms){
	enrich_out@result <- enrich_out@result[enrich_out@result$ID %in% terms,]
		
		return(enrich_out)
}

test2 <- enrichments_2p_h |> 
	map(keep, )
	map(filter_enrichment_by_term, common_terms_table) |> 
	plot_enrichment()

map(enrichments_2p_h, ~{.x[.x@result$ID]})


h_plots <- 
	h_enrichments |> 
	map(plot_enrichment) |> 
	imap(~{.x + labs(title = sample_id, subtitle = .y, caption = "h")})