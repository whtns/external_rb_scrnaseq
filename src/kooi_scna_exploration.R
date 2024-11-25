library(tidyverse)
library(plyranges)

kooi_scnas <- 
readxl::read_excel("data/kooi_somatic_genomic_2016/supp_4.xls", skip = 6) |> 
	dplyr::rename(seqnames = chrom, start= loc.start, end = loc.end) |> 
	dplyr::filter(seg.mean >= 1 | seg.mean <= -0.5) |>
	as_granges() |>
	identity()

arms_df <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
										col_names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")) |> 
	mutate(arm = substring(name, 1, 1)) |>
	group_by(chrom, arm) |>
	summarise(
		start = min(chromStart),
		end = max(chromEnd),
		length = end - start
	) |>
	dplyr::mutate(chrom = str_remove(chrom, "chr")) |>
	dplyr::mutate(chrom = dplyr::case_when(chrom == "X" ~ "23",
																				 chrom == "Y" ~ "24",
																				 .default = chrom)) |> 
	dplyr::filter(chrom %in% seq(1:24)) |> 
	dplyr::mutate(seqnames = chrom)

kooi_scna_ranges <-
	arms_df |>
	as_granges() |>
	join_overlap_intersect(kooi_scnas) |>
	as_tibble() |>
	dplyr::mutate(seqnames = str_pad(seqnames, 2, pad = "0")) |>
	dplyr::arrange(seqnames, arm) |>
	dplyr::select(seqnames, arm, everything()) |> 
	dplyr::rowwise() |> 
	dplyr::mutate(arm_percent = width/length) |> 
	dplyr::mutate(arm_effect = ifelse(arm_percent > 0.75, 1, 0)) |> 
	dplyr::group_by(ID) |>
	dplyr::summarise(n_events = sum(arm_effect)) |>
	dplyr::arrange(n_events) |>
	identity()
