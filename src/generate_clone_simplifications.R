#!/usr/bin/env Rscript
# Generate large_clone_simplifications YAML from numbat output
# Usage: Rscript src/generate_clone_simplifications.R SRR14800534

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

#' Parse genotype string and extract chromosome arms
#'
#' Converts genotypes like "1q+", "2p-", "13cnloh" to category labels
#' e.g. "1q", "2p", "13"
#'
#' @param genotypes Character vector of genotype strings (e.g. c("1q+", "2p-", "13cnloh"))
#' @return Data frame with columns: input_genotype, chrom, arm, category
parse_genotypes <- function(genotypes) {
  df <- tibble(input = genotypes) %>%
    # Extract chromosome number and arm/suffix
    mutate(
      # Match patterns like: "1q+", "1q-", "13cnloh", "21+", "5qcnloh"
      chrom = str_extract(input, "^[0-9]+"),
      arm = str_extract(input, "[pq](?=[+-]|cnloh|$)"),
      suffix_type = str_extract(input, "cnloh|[+-]|(?<=[pq])$"),
      # Build category: "1q", "2p", "13", etc.
      category = case_when(
        is.na(arm) & !is.na(chrom) ~ chrom,  # e.g. "21+" → "21"
        !is.na(arm) ~ paste0(chrom, arm),    # e.g. "1q+" → "1q"
        TRUE ~ NA_character_
      )
    ) %>%
    select(input_genotype = input, chrom, arm, suffix_type, category)
  
  return(df)
}

#' Generate YAML genotype→simplification mappings for a single sample
#'
#' @param sample_id SRR sample ID (e.g. "SRR14800534")
#' @param numbat_rds_path Path to the numbat RDS file
#' @return List (in YAML format) mapping genotypes to simplifications
#'
#' @details
#' For each unique genotype in a sample, extracts the chromosome arm(s) affected
#' and assigns a suffix letter (a, b, c, ...) based on clone ordering.
#'
#' Example output:
#'   SRR14800534:
#'     1q+: 1b
#'     16q-: 16c
#'
generate_simplifications_for_sample <- function(sample_id, numbat_rds_path) {
  cat(sprintf("\n=== Processing %s ===\n", sample_id))
  
  # Read numbat RDS
  if (!file.exists(numbat_rds_path)) {
    stop(sprintf("Numbat RDS not found: %s", numbat_rds_path))
  }
  
  mynb <- readRDS(numbat_rds_path)
  
  # Extract clone post-processing table
  clone_post <- mynb[["clone_post"]] %>%
    as_tibble(rownames = "cell") %>%
    select(cell, clone_opt, GT_opt) %>%
    # GT_opt may contain multiple genotypes separated by comma
    tidyr::separate_rows(GT_opt, sep = ",") %>%
    mutate(GT_opt = str_trim(GT_opt))
  
  cat(sprintf("Found %d cells, %d unique genotypes\n", 
              n_distinct(clone_post$cell),
              n_distinct(clone_post$GT_opt)))
  
  # Parse genotypes to extract arms
  genotypes_unique <- clone_post %>%
    select(GT_opt) %>%
    distinct() %>%
    pull(GT_opt)
  
  genotype_parsed <- parse_genotypes(genotypes_unique)
  
  cat("\nGenotype parsing:\n")
  print(genotype_parsed)
  
  # Map each genotype to (category, clone_opt) pairs
  # Order clones by first appearance in data (proxy for phylogenetic order)
  genotype_clone_map <- clone_post %>%
    select(clone_opt, GT_opt) %>%
    distinct() %>%
    left_join(genotype_parsed, by = c("GT_opt" = "input_genotype")) %>%
    # For each category, assign suffix letters based on clone order
    group_by(category) %>%
    arrange(clone_opt) %>%
    mutate(
      suffix = letters[seq_along(unique(clone_opt))],
      simplified = paste0(category, suffix)
    ) %>%
    ungroup()
  
  cat("\nGenotype→Clone→Simplified mapping:\n")
  print(genotype_clone_map)
  
  # Create genotype → simplification mapping
  # (Use first assigned simplification for each genotype)
  mapping <- genotype_clone_map %>%
    select(GT_opt, simplified) %>%
    distinct(GT_opt, .keep_all = TRUE) %>%
    deframe()  # Convert to named vector for YAML
  
  return(mapping)
}

# ─────────────────────────────────────────────────────────────────────────────

# MAIN
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript src/generate_clone_simplifications.R <sample_id>\n")
  cat("Example: Rscript src/generate_clone_simplifications.R SRR14800534\n")
  quit(status = 1)
}

sample_id <- args[1]

# Get the project root directory (parent of src/)
project_root <- dirname(dirname(getwd()))
if (!grepl("external_rb_scrnaseq_proj$", getwd())) {
  # Script called from outside project dir, use explicit path
  project_root <- "/home/skevin/single_cell_projects/resources/external_rb_scrnaseq_proj"
}

numbat_rds_path <- file.path(project_root, "output/numbat_sridhar", 
                              sprintf("%s_numbat.rds", sample_id))

# Generate mapping
mapping <- generate_simplifications_for_sample(sample_id, numbat_rds_path)

# Pretty print as YAML
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("YAML OUTPUT FOR large_clone_simplifications.yaml:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

yaml_output <- list(mapping)
names(yaml_output) <- sample_id

cat(as.yaml(yaml_output, indent.mapping = 2))
