#!/usr/bin/env Rscript

write_llm_context <- function(output_file = "targets_llm_context.txt", force = FALSE) {
  message("Gathering pipeline context...")
  
  # 1. Simple caching check: don't regenerate if metadata is newer than last export
  meta_path <- "_targets_r431/meta/meta"
  if (!force && file.exists(output_file) && file.exists(meta_path)) {
    if (file.info(meta_path)$mtime < file.info(output_file)$mtime) {
      message("Context file is up to date. Use force = TRUE to overwrite.")
      return(invisible(output_file))
    }
  }

  con <- file(output_file, "w")
  
  # --- Section 1: System & Pipeline Summary ---
  writeLines("=== PIPELINE SYSTEM SUMMARY ===", con)
  writeLines(paste("Export Date:", Sys.time()), con)
  writeLines(paste("Project Root:", getwd()), con)
  writeLines("\n", con)

  # --- Section 2: DAG Structure (Mermaid) ---
  writeLines("=== MERMAID DAG STRUCTURE ===", con)
  writeLines("```mermaid", con)
  writeLines(targets::tar_mermaid(), con)
  writeLines("```\n", con)

  # --- Section 3: Performance Metadata ---
  # Helps the LLM identify I/O bottlenecks on your JBOD
  writeLines("=== TARGETS METADATA (Performance & I/O) ===", con)
  if (dir.exists("_targets")) {
    meta <- targets::tar_meta()
    # Focus on timing and size to identify JBOD bottlenecks
    meta_clean <- meta[, intersect(names(meta), c("name", "type", "seconds", "bytes", "format"))]
    write.table(meta_clean, con, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  writeLines("\n", con)

  # --- Section 4: Code Definitions ---
  writeLines("=== CODE DEFINITIONS (_targets.R and R/ directory) ===", con)
  
  # Add _targets.R
  if (file.exists("_targets.R")) {
    writeLines("--- _targets.R ---", con)
    writeLines(readLines("_targets.R"), con)
  }
  
  # Add all functions in R/
  r_files <- list.files("R", full.names = TRUE, pattern = "\\.[Rr]$")
  for (f in r_files) {
    writeLines(paste("\n---", f, "---"), con)
    writeLines(readLines(f), con)
  }

  close(con)
  message(paste("Success! Context saved to:", output_file))
}

write_llm_context()