library(targets)
suppressPackageStartupMessages(source("packages.R"))
source("functions.R")
library(plotgardener)

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
	dplyr::arrange(desc(n_events)) |>
	identity()

# check kooi candidates ------------------------------

tar_load("table_all_diffex_clones")

col_types <- c("guess", "guess", "guess", "guess", "guess", "guess", "guess", 
							 "guess", "guess", "guess", "guess", "guess", "guess", "guess", 
							 "guess", "guess", "guess", "guess", "guess", "guess", "guess", 
							 "guess", "text", "numeric")

total_kooi_candidates <- 
	c("1", "6", "16") |> 
	set_names() |> 
	map(\(x) read_excel(path = "results/diffex_bw_clones_all_by_chr.xlsx", sheet = x, col_types = col_types)) |> 
	map(dplyr::filter, !is.na(kooi_region)) |> 
	dplyr::bind_rows() |> 
	dplyr::distinct(sample_id, symbol, p_val, .keep_all = TRUE) |> 
	dplyr::mutate(kept = dplyr::case_when(
		(chr == "1" & avg_log2FC > 0) ~ 1,
		(chr == "6" & avg_log2FC > 0) ~ 1,
		(chr == "16" & avg_log2FC < 0) ~ 1,
		.default = 0
	)) |> 
	# dplyr::filter(kept == 1) |> 
	dplyr::arrange(location) |> 
identity()


col_types <- c("guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", 
							 "guess", "guess", "guess", "guess", "guess", "guess", "guess", 
							 "guess", "guess", "guess", "guess", "guess", "guess", "guess", 
							 "text", "text", "numeric")

cluster_kooi_candidates <- 
	c("1", "6", "16") |> 
	set_names() |> 
	map(\(x) read_excel(path = "results/diffex_bw_clones_per_cluster_all_by_chr.xlsx", sheet = x, col_types = col_types)) |>
	map(dplyr::filter, !is.na(kooi_region)) |> 
	dplyr::bind_rows() |>
	dplyr::distinct(sample_id, symbol, p_val, .keep_all = TRUE) |> 
	dplyr::mutate(kept = dplyr::case_when(
		chr == "1" & avg_log2FC > 0 ~ 1,
		chr == "6" & avg_log2FC > 0 ~ 1,
		chr == "16" & avg_log2FC < 0 ~ 1,
		.default = 0
	)) |> 
	# dplyr::filter(kept == 1) |> 
	identity()

writexl::write_xlsx(list("whole tumor" = total_kooi_candidates, "cluster" = cluster_kooi_candidates), "results/table_s02a.xlsx") |> 
	browseURL()


# mean expression of putative drivers ------------------------------

kooi_candidates <- read_csv("data/kooi_candidates.csv")

tar_load("final_seus")

pull_gene_expression_from_seu <- function(seu_path){
	seu <- readRDS(seu_path)
	
	features <- kooi_candidates[["symbol"]]
	
	expr_mat <- GetAssayData(seu, assay = "gene")
	
	rowMeans(expr_mat[rownames(expr_mat) %in% kooi_candidates[["symbol"]],]) |> 
		tibble::enframe("symbol", "expression")
	
}

kooi_candidate_mean_expression <- map(final_seus, pull_gene_expression_from_seu)

test0 <- bind_rows(kooi_candidate_mean_expression, .id = "sample_id") |> 
	dplyr::left_join(kooi_candidates, by = "symbol") |> 
	dplyr::mutate(kooi_region = factor(kooi_region, levels = c("1q", "2p", "6p", "16q")))

write_csv(test0, "results/kooi_candidate_mean_expression.csv")

browseURL("results/kooi_candidate_mean_expression.csv")

test0 |> 
	dplyr::mutate(symbol = factor(symbol, levels = unique(kooi_candidates$symbol))) |> 
	group_by(symbol, kooi_region) |> 
	dplyr::summarise(mean_expression = mean(expression)) |> 
	# View() |> 
	dplyr::pull(mean_expression) |>
	cat(sep = "\n") |>
	identity()
	


ggplot(test0, aes(x = symbol, y = expression, fill = kooi_region)) +
	geom_boxplot() + 
	facet_wrap(~kooi_region, scales = "free") + 
	ylim(0, 4) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
	labs(fill = "SCNA") + 
	NULL

ggsave("results/kooi_candidate_mean_expression.pdf", w = 8, h = 4) |> 
	browseURL()


