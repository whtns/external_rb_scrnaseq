## Load Giemsa stain band information and genomic
## annotation data for hg19 genome assembly
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationHub)
library(plotgardener)
library(readr)
library(stringr)
library(glue)
library(fs)
library(tidyverse)

cytobands <- suppressMessages(AnnotationHub::query(AnnotationHub::AnnotationHub(), 
																									 "AHCytoBands"))
cytobands

cytoData <-
	as.data.frame(suppressMessages(cytobands[["AH53178"]]))

make_annoHighlight_from_consensus <- function(ideogramPlot, ploty, chrom, chromstart, chromend, fill){
	region <- pgParams(chrom = chrom, chromstart = chromstart, chromend = chromend)
	
	annoHighlight(
		plot = ideogramPlot, params = region,
		fill = fill,
		y = ploty-0.125, height = 0.25, just = c("left", "top"), default.units = "inches"
	)
	
	if(chrom == "chr2"){
		chrom_lengths = seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)
		genesPlot <- plotGenes(
			chrom = "chr2",
			chromstart = 1,
			chromend = chrom_lengths[mychrom],
			assembly = "hg38",
			width = 6.25,
			geneHighlights = data.frame(
				"gene" = c("MYCN"),
				"color" = c("#225EA8")
			),
			geneBackground = NA,
			x = 0.25, y = 2.25, height = 0.75,
			# just = c("left", "top"), 
			default.units = "inches"
		)
	}
	
}

make_rb_scna_ideograms <- function(consensus_file) {
	
	tumor_id <- str_extract(consensus_file, "SRR[0-9]*")
	
	plot_path <- glue("results/{tumor_id}_karyogram.pdf")
	
  test0 <- read_tsv(consensus_file) |> 
  	dplyr::mutate(CHROM= paste0("chr", CHROM)) |> 
  	dplyr::filter(!is.na(seg)) |> 
  	dplyr::mutate(fill = dplyr::case_when(cnv_state == "amp" ~ "red",
  																				cnv_state == "del" ~ "blue")) |> 
  	dplyr::filter(fill %in% c("blue", "red")) |> 
  	dplyr::rename(chrom = CHROM, chromstart = seg_start, chromend = seg_end) |> 
  	dplyr::select(chrom, chromstart, chromend, fill)
  
  
  pdf(plot_path)
  
  pageCreate(
  	width = 6.25, height = 5.25, default.units = "inches",
  	showGuides = FALSE, xgrid = 0, ygrid = 0
  )
  plotText(
  	label = tumor_id, 
  	x = 1, y = 0.5,
  	fontsize = 20
  )
  
  mychrom = "chr1"
  ploty = 1.5
  # chr1------------------------------
  ideogramPlot <- plotIdeogram(
  	chrom = mychrom, assembly = "hg38",
  	orientation = "h",
  	x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
  )
  plotText(
  	label = gsub("chr", "", mychrom),
  	x = 0.25, y = ploty - 0.25, fontsize = 10
  )
  
  if(mychrom %in% test0$chrom){
  	print("yes")
  	test0 |> 
  		dplyr::filter(chrom %in% mychrom) |> 
  		pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
  }
  
  
  # chr2 ------------------------------
  mychrom = "chr2"
  ploty = 2.5
  ideogramPlot <- plotIdeogram(
  	chrom = mychrom, assembly = "hg38",
  	orientation = "h",
  	x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
  )
  plotText(
  	label = gsub("chr", "", mychrom),
  	x = 0.25, y = ploty - 0.25, fontsize = 10
  )
  
  if(mychrom %in% test0$chrom){
  	print("yes")
  	test0 |>
  		dplyr::filter(chrom %in% mychrom) |>
  		pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
  }
  

  # chr6 ------------------------------
  mychrom = "chr6"
  ploty = 3.5
  ideogramPlot <- plotIdeogram(
  	chrom = mychrom, assembly = "hg38",
  	orientation = "h",
  	x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
  )
  plotText(
  	label = gsub("chr", "", mychrom),
  	x = 0.25, y = ploty - 0.25, fontsize = 10
  )
  
  if(mychrom %in% test0$chrom){
  	print("yes")
  	test0 |> 
  		dplyr::filter(chrom %in% mychrom) |> 
  		pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
  }
  
  # chr16 ------------------------------
  mychrom = "chr16"
  ploty = 4.5
  ideogramPlot <- plotIdeogram(
  	chrom = mychrom, assembly = "hg38",
  	orientation = "h",
  	x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
  )
  plotText(
  	label = gsub("chr", "", mychrom),
  	x = 0.25, y = ploty - 0.25, fontsize = 10
  )
  
  if(mychrom %in% test0$chrom){
  	print("yes")
  	test0 |> 
  		dplyr::filter(chrom %in% mychrom) |> 
  		pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
  }
  
  dev.off()
  
  return(plot_path)
}

consensus_segs <- dir_ls("output/numbat_sridhar/", regexp = ".*[0-9]\\/segs_consensus_2.tsv", recurse = TRUE) |> 
	sort()

names(consensus_segs) <- str_extract(consensus_segs, "SRR[0-9]*")

make_rb_scna_ideograms(consensus_segs[[3]]) |> 
	browseURL()

plot_fig_s0x()

consensus_segs[c("SRR13884246", "SRR13884247", "SRR13884248", "SRR17960484")] |> 
	map(make_rb_scna_ideograms) |> 
	map(browseURL)
	

