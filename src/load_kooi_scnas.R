library(data.table)

chrom_coords <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz",
           col.names = c("chrom","chromStart","chromEnd","name","gieStain")) %>%
  dplyr::mutate(arm = dplyr::case_when(
    str_detect(name, "p") ~ "p",
    str_detect(name, "q") ~ "q"
  )) %>%
  dplyr::group_by(chrom, arm) %>%
  dplyr::summarize(start = min(chromStart), end = max(chromEnd)) %>%
  dplyr::mutate(seqnames = str_remove(chrom, "chr")) %>%
  dplyr::select(-chrom) %>%
  plyranges::as_granges() %>%
  identity()

kooi_gene_measures <- "data/kooi_scna_meta_analysis/kooi_dna.csv" %>%
  read_csv() %>%
  identity()

colnames(kooi_gene_measures) <- c("geneSymbol", 1:56)

kooi_scnas <-
  annotables::grch38 %>%
  dplyr::left_join(kooi_gene_measures, by = c("symbol" = "geneSymbol")) %>%
  dplyr::filter(chr %in% c("1", "2", "6", "16")) %>%
  dplyr::arrange(chr, start) %>%
  tidyr::pivot_longer(-colnames(annotables::grch38), names_to = "sample_id", values_to = "copyNumber") %>%
  drop_na()

ggplot(kooi_scnas, aes(start, copyNumber, color = sample_id)) +
  geom_line() +
  facet_wrap(~chr)

test0 <-
  kooi_scnas %>%
  dplyr::mutate(seqnames = chr) %>%
  plyranges::as_granges() %>%
  plyranges::join_overlap_intersect(chrom_coords) %>%
  as_tibble() %>%
  # dplyr::left_join(chrom_coords, by = "chr", relationship = "many-to-many") %>%
  tidyr::unite(chr_arm, seqnames, arm) %>%
  group_by(chr_arm, sample_id) %>%
  dplyr::summarize(mean_cn = mean(copyNumber)) %>%
  dplyr::filter(!between(mean_cn, -0.2, 0.05)) %>%
  identity()

ggplot(test0, aes(chr_arm, mean_cn, color = chr_arm)) +
  geom_point() +
  geom_hline(aes(yintercept = 0.05))

stachelek_scores = read_csv("results/stachelek_scores.csv") %>%
  tidyr::separate(name, c("region", "chr_arm"), "_") %>%
  dplyr::mutate(chr_arm = str_remove(chr_arm, "\\+|\\-")) %>%
  dplyr::mutate(chr_arm = str_replace(chr_arm, "q", "_q")) %>%
  dplyr::mutate(chr_arm = str_replace(chr_arm, "p", "_p")) %>%
  dplyr::rename(geneSymbol = value) %>%
  # split(.$name) %>%
  # purrr::map(pull, value) %>%
  identity()

scna_samples <-
  test0 %>%
  dplyr::filter(chr_arm %in% unique(stachelek_scores$chr_arm)) %>%
  split(., .$chr_arm) %>%
  # dplyr::bind_rows(.id = "region") %>%
  # purrr::map(dplyr::pull, sample_id) %>%
  identity()

kooi_rna_input <- read_csv("data/kooi_scna_meta_analysis/kooi_rna.csv")
colnames(kooi_rna_input) <- c("geneSymbol", "rnaProbeId", 1:56)
kooi_rna <-
  kooi_rna_input %>%
  dplyr::left_join(annotables::grch38, by = c("geneSymbol" = "symbol")) %>%
  tidyr::pivot_longer(-any_of(c("geneSymbol", "rnaProbeId", colnames(annotables::grch38))), names_to = "sample_id", values_to = "expression") %>%
  identity()

pdf("results/bulk_segment_scores.pdf")
for(segment in names(scna_samples)){

  stachelek_genes <-
    stachelek_scores %>%
    dplyr::filter(chr_arm == segment)

  res <-
    kooi_rna %>%
    dplyr::left_join(stachelek_genes, by = "geneSymbol") %>%
    dplyr::left_join(scna_samples[[segment]], by = c("sample_id", "chr_arm")) %>%
    dplyr::group_by(sample_id, chr_arm, region) %>%
    dplyr::summarize(signature = mean(expression)) %>%
    identity()

  myplot <- ggplot(res, aes(chr_arm, signature)) +
    facet_wrap(~region) +
    geom_boxplot() +
    labs(title = segment)

  print(myplot)

}
dev.off()

browseURL("results/bulk_segment_scores.pdf")
