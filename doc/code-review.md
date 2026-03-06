Files Reviewed: seu_metadata_db.R, integration_functions.R, plot_functions_12.R, metadata_functions_1.R

CRITICAL Issues
C1 — eval(parse(text = ...)) — code injection risk

integration_functions.R:59: with(seu@meta.data, eval(parse(text = filter_expr))) — caller-controlled string executed as R code
metadata_functions_1.R:187: eval(parse(text = x)) fallback when JSON parsing fails on a value read from SQLite — stored-code-injection vector
Fix: Replace filter_expr string eval with a safe typed mechanism (e.g. named list of column/value conditions). Remove the eval fallback in read_cluster_orders_table() entirely — migrate any legacy non-JSON rows.

HIGH Issues
H1 — Dead code / broken contract — metadata_functions_1.R:47-49, 73-75
Both add_batch_hash_metadata() and add_hash_metadata() have an unreachable return(seu) after return(filepath). Roxygen says they return a Seurat object — they actually return a filepath string.

H2 — dbDisconnect() not in on.exit() — metadata_functions_1.R
seu_metadata_db.R uses on.exit(DBI::dbDisconnect(con)) correctly everywhere. metadata_functions_1.R does not — connections leak on any error or early return(). In read_batch_hashes(), the dbDisconnect() is placed after an early return() and is never reached on that code path.

H3 — compare_cluster_composition() uses glue() for SQL column interpolation — seu_metadata_db.R:342-351
group_by is guarded by match.arg() now, but glue interpolation into SQL is inherently unsafe. Use switch() to select a pre-written query string instead.

H4 — No tests for any reviewed functions
Tests exist only for transcripts_to_genes. None of the DB functions, metadata functions, or plot helpers have coverage.

H5 — retrieve_snakemake_params() references undefined sample_id — metadata_functions_1.R:220
str_extract() result is discarded (not assigned), then sample_id is used in return(list(sample_id, params)) — this throws Error: object 'sample_id' not found at runtime.

H6 — seu_metadata_db.R is 401 lines — marginally over the 400-line threshold; mixes schema init, extraction, and query helpers.

MEDIUM Issues
M1 — Pervasive %>% instead of |> — 2,419 occurrences across 62 files. Inconsistently mixed with native pipe in newer code.

M2 — select_genes_to_plot() has no explicit return() — plot_functions_12.R:42-106. Last-expression return differs between branches — fragile.

M3 — Commented-out browser() calls — integration_functions.R:44 and several other files.

M4 — Old-style join syntax by = c("a" = "b") — should use join_by(a == b) in multiple files.

M5 — group_by() + summarize() without .by — plot_functions_12.R:14-15.

M6 — make_filepaths_unique_in_hashes_table() may not collect before overwrite — metadata_functions_1.R:93-110. Needs explicit collect() + transaction.

M7 — Hardcoded relative paths ("output/seurat/...") in multiple functions — should use here::here().

LOW Issues
L1 — Placeholder roxygen descriptions in plot_functions_12.R:32-41.
L2 — integration_by_scna_clones() has no documentation.
L3 — possibly_plot_clone_cc_plots assigned twice identically in plot_functions_12.R:147,171.
L4 — sapply() used in 10+ places (type-unstable; use map_*()).
L5 — DEFAULT_DB is a package-level global at seu_metadata_db.R:12.

Verdict: NEEDS CHANGES
Severity	Count
CRITICAL	2
HIGH	6
MEDIUM	7
LOW	5
Required before committing: Fix C1 (eval/parse injection), H1 (dead code/broken return contract), H2 (dbDisconnect leak), H5 (undefined sample_id runtime error).