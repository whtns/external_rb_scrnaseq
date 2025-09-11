#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

genes_tbl <- annotables_to_grange()

depmap_rnai_by_range <- function(){
	
	arms_granges <- get_arms_ranges()
	
	rb_scna_granges <- arms_granges[arms_granges$chrom_arm %in% c("1q", "2p", "6p", "16q")]
	
	depmap_rnai_granges <- 
		depmap::depmap_rnai() |> 
		# dplyr::mutate(symbol = str_extract("A1BG (1)", ".*(?=\\s)")) |> 
		dplyr::left_join(annotables::grch38, by = c("gene_name" = "symbol")) |> 
		dplyr::mutate(seqnames = chr) |> 
		plyranges::as_granges() |> 
		identity()
	
}

drna <- depmap::depmap_rnai()

