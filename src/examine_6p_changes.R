library(targets)
suppressPackageStartupMessages(source("packages.R"))
source("functions.R")
library(plotgardener)
library(plyranges)

tar_load(c("corresponding_clusters_diffex_6p", "corresponding_state_6p_seus", "corresponding_clusters_enrichments_6p"))

segs_6p <- sapply(c("output/numbat_sridhar/SRR13884247_numbat.rds",
										"output/numbat_sridhar/SRR17960484_numbat.rds"), pull_scna_segments, chrom = "6") |> 
	as("GRangesList") |> 
	unlist() |> 
	# reduce_ranges() |> 
	identity()



test2 <- corresponding_clusters_enrichments_6p$corresponding_clusters_enrichments_6p_0c262e28$enrichment$`g1_0 v. g1_1` |> 
	clusterProfiler::setReadable(OrgDb = 'org.Hs.eg.db', keyType = "ENTREZID")

test3 <- 
	test2@result |> 
	dplyr::filter(pvalue <= 0.1) |> 
	dplyr::select(Description, core_enrichment) %>%
	dplyr::mutate(core_enrichment = str_split(core_enrichment, pattern = "/")) %>%
	tidyr::unnest(core_enrichment) %>%
	dplyr::left_join(annotables::grch38, by = c("core_enrichment" = "symbol")) %>%
	dplyr::mutate(seqnames = chr) |> 
	dplyr::filter(!is.na(start), !is.na(end)) |> 
	as_granges() |> 
	identity()

test4 <-
	test3 |> 
	setdiff_ranges(segs_6p[1,]) |> 
	join_overlap_intersect(test3) |> 
	as_tibble() |>
	dplyr::filter(seqnames == "6") |> 
	# dplyr::filter(avg_log2FC > 0) |> 
	identity()

names(corresponding_clusters_diffex_6p) <- str_extract(corresponding_state_6p_seus, "SRR[0-9]*")

cis_47 <- corresponding_clusters_diffex_6p$SRR13884247$`g1_0-g1_1 v. g1_4-g1_5`$cis

cis_84 <- corresponding_clusters_diffex_6p$SRR17960484$`g1_0 v. g1_1`$cis

specific_47 <- anti_join(cis_47, cis_84, by = "symbol") |> 
	as_granges() |> 
	setdiff_ranges(segs_6p[2,])

specific_84 <- 
	cis_84 |> 
	# anti_join(cis_47, by = "symbol") |> 
	as_granges() |> 
	identity()

final_84 <-
	specific_84 |> 
	setdiff_ranges(segs_6p[1,]) |> 
	join_overlap_intersect(specific_84) |> 
		as_tibble() |> 
	dplyr::filter(avg_log2FC > 0) |> 
		identity()

View(enrichment_analysis)

final_84 |> 
	tibble::column_to_rownames("symbol") %>%
	# dplyr::select(-any_of(c("entrez"))) |> 
	# enrichment_analysis(analysis_method = "ora", pvalueCutoff = 1) |>
	dplyr::pull(entrez) |> 
	clusterProfiler::enrichGO(
		OrgDb = 'org.Hs.eg.db',
		ont = "BP",
		pvalueCutoff = 1,
		pAdjustMethod = "BH"
	) |> 
	clusterProfiler::dotplot() |>
	identity()

# test0 ------------------------------

gene_sets <- msigdbr::msigdbr(species = "human", category = "H") |> 
	dplyr::select(gs_name, entrez_gene)

gse <- clusterProfiler::enricher(
	gene = gene_list,
	minGSSize = 3,
	maxGSSize = 800,
	pvalueCutoff = 1,
	TERM2GENE = gene_sets,
	pAdjustMethod = "BH"
)

# test0 ------------------------------


library(msigdbr)

gene_set_genes <- msigdbr::msigdbr(category = "H") |> 
	dplyr::filter(gs_name == "HALLMARK_MYC_TARGETS_V1") |> 
	dplyr::pull(gene_symbol)

final_84 |> dplyr::filter(symbol %in% gene_set_genes)
