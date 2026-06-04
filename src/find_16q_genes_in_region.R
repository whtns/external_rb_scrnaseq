library(tidyverse)

# centromere ------------------------------
centromere_16q_genes <-
  annotables::grch38 %>%
  dplyr::arrange(chr, start) %>%
  dplyr::filter(chr == "16") %>%
  dplyr::filter(start >= 46469341) %>%
  dplyr::filter(!str_detect(symbol, "^AC0")) %>%
  dplyr::pull(symbol) %>%
  unique() %>%
  head(10) %>%
  cat(sep = "\n")

centro_grnas <- read_tsv("data/gRNAs/centro_gRNAs_raw.txt") %>%
  janitor::clean_names() %>%
  dplyr::left_join(annotables::grch38, by = c("target_gene_symbol" = "symbol")) %>%
  dplyr::mutate(bed_start = sg_rna_cut_position_1_based) %>%
  dplyr::mutate(bed_end = bed_start + nchar(sg_rna_sequence)) %>%
  dplyr::mutate(sign = dplyr::case_when(
    strand > 0 ~ "+",
    strand < 0 ~ "-")) %>%
  dplyr::mutate(strand = "0") %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  dplyr::select(chr, bed_start, bed_end, name, strand, sign, everything()) %>%
  dplyr::filter(target_cut_percent > 70) %>%
  dplyr::arrange(bed_start) %>%
  dplyr::mutate(grna_diff = dplyr::lead(bed_start) - bed_start) %>%
  dplyr::mutate(lead_name = dplyr::lead(name)) %>%
  dplyr::select(grna_diff, lead_name, everything()) %>%
  dplyr::arrange(grna_diff) %>%
  identity()

centro_grnas %>%
  dplyr::select(chr, bed_start, bed_end, name, strand, sign, everything()) %>%
  write_delim("data/gRNAs/centro_grnas.bed", col_names = FALSE)

# telomere ------------------------------
telomere_16q_genes <-
  annotables::grch38 %>%
  dplyr::arrange(chr, desc(start)) %>%
  dplyr::filter(chr == "16") %>%
  dplyr::filter(start >= 89000000) %>%
  dplyr::filter(!str_detect(symbol, "^AC")) %>%
  dplyr::filter(!str_detect(symbol, "^MIR")) %>%
  dplyr::filter(!str_detect(symbol, "^RNU")) %>%
  dplyr::filter(!str_detect(symbol, "^LINC")) %>%
  dplyr::pull(symbol) %>%
  unique() %>%
  head(10) %>%
  cat(sep = "\n")

telo_grnas <- read_tsv("data/gRNAs/telo_gRNAs_raw.txt") %>%
  janitor::clean_names() %>%
  dplyr::left_join(annotables::grch38, by = c("target_gene_symbol" = "symbol")) %>%
  dplyr::mutate(bed_start = sg_rna_cut_position_1_based) %>%
  dplyr::mutate(bed_end = bed_start + nchar(sg_rna_sequence)) %>%
  dplyr::mutate(sign = dplyr::case_when(
    strand > 0 ~ "+",
    strand < 0 ~ "-")) %>%
  dplyr::mutate(strand = "0") %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  dplyr::select(chr, bed_start, bed_end, name, strand, sign, everything()) %>%
  dplyr::filter(target_cut_percent > 70) %>%
  dplyr::arrange(bed_start) %>%
  dplyr::mutate(grna_diff = dplyr::lead(bed_start) - bed_start) %>%
  dplyr::mutate(lead_name = dplyr::lead(name)) %>%
  dplyr::select(grna_diff, lead_name, everything()) %>%
  dplyr::arrange(grna_diff) %>%
  identity()

telo_grnas %>%
  dplyr::select(chr, bed_start, bed_end, name, strand, sign, everything()) %>%
  write_delim("data/gRNAs/telo_grnas.bed", col_names = FALSE)


# midpoint ------------------------------
midpoint_16q = 46469341+(90000000-46469341)/2

midpoint_16q_genes <-
  annotables::grch38 %>%
  dplyr::arrange(chr, desc(start)) %>%
  dplyr::filter(chr == "16") %>%
  dplyr::filter(between(start, midpoint_16q-5e5, midpoint_16q+5e5)) %>%
  dplyr::filter(!str_detect(symbol, "^AC")) %>%
  dplyr::filter(!str_detect(symbol, "^MIR")) %>%
  dplyr::filter(!str_detect(symbol, "^RNU")) %>%
  dplyr::pull(symbol) %>%
  unique() %>%
  head(20) %>%
  cat(sep = "\n")


midpoint_grnas <- read_tsv("data/gRNAs/midpoint_gRNAs_raw.txt") %>%
  janitor::clean_names() %>%
  dplyr::left_join(annotables::grch38, by = c("target_gene_symbol" = "symbol")) %>%
  dplyr::mutate(bed_start = sg_rna_cut_position_1_based) %>%
  dplyr::mutate(bed_end = bed_start + nchar(sg_rna_sequence)) %>%
  dplyr::mutate(sign = dplyr::case_when(
    strand > 0 ~ "+",
    strand < 0 ~ "-")) %>%
  dplyr::mutate(strand = "0") %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  dplyr::select(chr, bed_start, bed_end, name, strand, sign, everything()) %>%
  dplyr::filter(target_cut_percent > 70) %>%
  dplyr::arrange(bed_start) %>%
  dplyr::mutate(grna_diff = dplyr::lead(bed_start) - bed_start) %>%
  dplyr::mutate(lead_name = dplyr::lead(name)) %>%
  dplyr::select(grna_diff, lead_name, everything()) %>%
  dplyr::arrange(grna_diff) %>%
  identity()

midpoint_grnas %>%
  dplyr::select(chr, bed_start, bed_end, name, strand, sign, everything()) %>%
  write_delim("data/gRNAs/midpoint_grnas.bed", col_names = FALSE)




