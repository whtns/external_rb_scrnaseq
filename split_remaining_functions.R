#!/usr/bin/env Rscript

# Script to complete the splitting of functions.R into multiple files
# This script extracts all remaining functions and splits them into files with max 5 functions each

library(stringr)
library(readr)

# Read the original functions.R file
functions_content <- read_file("R/functions.R")

# Extract all function definitions using regex
function_pattern <- "([a-zA-Z_][a-zA-Z0-9_.]*\\s*<-\\s*function\\s*\\([^)]*\\)\\s*\\{)"
function_matches <- str_locate_all(functions_content, function_pattern)[[1]]

if (nrow(function_matches) == 0) {
  cat("No functions found.\n")
  quit()
}

# Extract function names
function_names <- str_extract_all(functions_content, function_pattern)[[1]]
function_names <- str_extract(function_names, "^[a-zA-Z_][a-zA-Z0-9_.]*")

# Find the start and end of each function
function_starts <- function_matches[, "start"]
function_ends <- numeric(length(function_starts))

# For each function, find where it ends by counting braces
for (i in seq_along(function_starts)) {
  start_pos <- function_starts[i]
  brace_count <- 0
  pos <- start_pos
  in_function_def <- FALSE
  
  while (pos <= nchar(functions_content)) {
    char <- substr(functions_content, pos, pos)
    
    if (char == "{") {
      brace_count <- brace_count + 1
      in_function_def <- TRUE
    } else if (char == "}") {
      brace_count <- brace_count - 1
      if (in_function_def && brace_count == 0) {
        function_ends[i] <- pos
        break
      }
    }
    pos <- pos + 1
  }
}

# Extract each complete function
functions_list <- list()
for (i in seq_along(function_starts)) {
  if (function_ends[i] > 0) {
    func_text <- substr(functions_content, function_starts[i], function_ends[i])
    functions_list[[function_names[i]]] <- func_text
  }
}

# Group functions into files of 5 each
functions_per_file <- 5
function_groups <- split(names(functions_list), ceiling(seq_along(functions_list) / functions_per_file))

# Create themed file names based on function content
get_file_theme <- function(func_names) {
  func_text <- paste(functions_list[func_names], collapse = " ")
  
  if (str_detect(func_text, "plot|ggplot|FeaturePlot|DimPlot")) {
    return("plot_functions")
  } else if (str_detect(func_text, "diffex|FindMarkers|enrichment")) {
    return("diffex_functions")
  } else if (str_detect(func_text, "numbat|clone|heatmap")) {
    return("numbat_functions")
  } else if (str_detect(func_text, "cluster|seurat|seu")) {
    return("cluster_functions")
  } else if (str_detect(func_text, "score|pca|qc|filter")) {
    return("scoring_functions")
  } else if (str_detect(func_text, "read|write|save|load")) {
    return("io_functions")
  } else {
    return("utility_functions")
  }
}

# Check existing files to avoid overwriting
existing_files <- list.files("R/", pattern = "\\.R$")
existing_count <- length(existing_files)

# Create files for remaining function groups
for (i in seq_along(function_groups)) {
  group_names <- function_groups[[i]]
  theme <- get_file_theme(group_names)
  
  # Find next available number for this theme
  theme_files <- str_subset(existing_files, paste0("^", theme))
  if (length(theme_files) > 0) {
    theme_numbers <- str_extract(theme_files, "\\d+") %>% as.numeric() %>% max(na.rm = TRUE)
    file_number <- theme_numbers + 1
  } else {
    file_number <- 1
  }
  
  filename <- paste0("R/", theme, "_", file_number, ".R")
  
  # Skip if we've already created files with similar names in this session
  if (file.exists(filename)) {
    filename <- paste0("R/", theme, "_", existing_count + i, ".R")
  }
  
  # Create file header
  file_content <- paste0("# ", str_to_title(str_replace_all(theme, "_", " ")), " (", existing_count + i, ")\n\n")
  
  # Add each function to the file
  for (func_name in group_names) {
    if (func_name %in% names(functions_list)) {
      file_content <- paste0(file_content, functions_list[[func_name]], "\n\n")
    }
  }
  
  # Write the file
  write_file(file_content, filename)
  cat("Created:", filename, "with functions:", paste(group_names, collapse = ", "), "\n")
}

cat("\nFunction splitting completed!\n")
cat("Total files created:", length(function_groups), "\n")
cat("Use source('R/load_all_functions.R') to load all functions.\n")