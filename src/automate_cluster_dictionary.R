#!/usr/bin/env Rscript
# Regenerate the automated cluster dictionary + discrepancy report using the
# same numbatHelpers::build_cluster_dictionary() that the targets pipeline now
# calls. This is the standalone/inspection entry point.
#
# Outputs (results/cluster_dictionary/):
#   cluster_dictionary_auto.tsv          (drop-in replacement candidate)
#   cluster_dictionary_discrepancies.tsv (old vs new, per cluster)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)

# base unfiltered seus, named by sample id (same set as the all_seu_files target)
seu_files <- fs::dir_ls("output/seurat/", regexp = "/SRX[0-9]+_seu.rds")
names(seu_files) <- str_extract(seu_files, "SRX[0-9]+")
# restrict to samples present in the existing hand dictionary, for a clean diff
dict_samples <- read_tsv("data/cluster_dictionary.tsv", show_col_types = FALSE)$sample_id |> unique()
seu_files <- seu_files[names(seu_files) %in% dict_samples]

dict <- build_cluster_dictionary(
  seu_files,
  sqlite_path     = "batch_hashes.sqlite",
  malat1_rank_max = 5L,
  top_n           = 10L,
  remove_programs = c("rod", "cone", "APOE", "low_qual", "MALAT1"),
  write_audit_dir = "results/cluster_dictionary",
  compare_to      = "data/cluster_dictionary.tsv"
)

rep <- read_tsv("results/cluster_dictionary/cluster_dictionary_discrepancies.tsv",
                show_col_types = FALSE)

cat("\n================ MALAT1 clusters: old -> new ================\n")
print(as.data.frame(rep |>
  filter(old_abbreviation == "MALAT1" | abbreviation == "MALAT1") |>
  select(sample_id, gene_snn_res.0.2, old_abbreviation, abbreviation, rule, remove)))

cat("\n================ new abbreviation distribution ================\n")
print(rep |> count(abbreviation, sort = TRUE), n = 60)

cat("\n================ clusters now flagged remove=1, by program ================\n")
print(rep |> filter(remove == "1") |> count(abbreviation, sort = TRUE))

cat("\nClusters changed: ", sum(rep$changed, na.rm = TRUE), " / ", nrow(rep), "\n")
cat("Samples annotated: ", length(dict), "\n")
