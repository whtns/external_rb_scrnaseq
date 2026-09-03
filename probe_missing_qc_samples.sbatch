#!/bin/bash
#SBATCH --job-name=probe_qc
#SBATCH --account=cobrinik_1090
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/probe_qc_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/probe_qc_%j.err

# Four SRX samples have NO rows in cell_qc_values, so they would get no depth
# panel (github #43): SRX10031193, SRX10264520, SRX10264521, SRX10264522.
# Their _seu.rds all exist (356 MB - 1.4 GB).
#
# bulk_extract_seu_metadata() is purrr::walk() with no error handling, so one
# bad object fails the whole seu_metadata_db target. Before spending 24 h on a
# production run, find out whether these four can actually be ingested.
#
# This DOES write to batch_hashes.sqlite (idempotent upsert, the same call the
# seu_metadata_db target makes). Nothing else may be running against the store.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)

missing <- c("SRX10031193", "SRX10264520", "SRX10264521", "SRX10264522")

for (s in missing) {
  f <- file.path("output/seurat", paste0(s, "_seu.rds"))
  cat("\n=== ", s, " ===\n", sep = "")
  t0 <- Sys.time()
  # one at a time, guarded -- the production target is not guarded, so a failure
  # here is exactly what would take down seu_metadata_db
  err <- tryCatch({
    extract_seu_metadata(f, sqlite_path = "batch_hashes.sqlite")
    NULL
  }, error = function(e) conditionMessage(e))
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (is.null(err)) {
    d <- get_sample_cell_depth(s)
    cat(sprintf("  ingested OK in %ds; depth rows now: %d", el, nrow(d)))
    if (nrow(d)) cat(sprintf("  median=%.0f min=%.0f", median(d$nCount_gene), min(d$nCount_gene)))
    cat("\n")
  } else {
    cat("  FAILED after ", el, "s: ", err, "\n", sep = "")
  }
}

cat("\n=== coverage after probe ===\n")
con <- connect_hash_db("batch_hashes.sqlite"); on.exit(DBI::dbDisconnect(con))
n <- DBI::dbGetQuery(con, "SELECT COUNT(DISTINCT sample_id) n FROM cell_qc_values WHERE sample_id LIKE \"SRX%\"")$n
cat("distinct SRX samples in cell_qc_values: ", n, " (cohort is 39)\n", sep = "")
'
