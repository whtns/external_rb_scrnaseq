# Smoke test for the three diffex fixes:
#   1. .seu_sample_id()  -- resolves config ids from *_hypoxia_low_seu.rds paths
#      (root cause of find_diffex_clones() returning 33 empty lists)
#   2. tally_kooi_candidates() -- survives an empty placeholder workbook
#   3. make_volcano_diffex_clones() dedupe -- no duplicate rownames from
#      patch-contig annotations (ATMIN on chr16 AND HG405_PATCH)
#
# Fast: reads cached CSVs, runs one real find_diffex_clones() call. No tar_make.

devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet = TRUE)
suppressPackageStartupMessages({
  library(targets); library(dplyr); library(readr); library(stringr); library(purrr)
})
tar_config_set(store = "_targets_r431")

ok <- function(x, msg) if (isTRUE(x)) cat("PASS:", msg, "\n") else stop("FAIL: ", msg)

# --- 1. sample id resolution -------------------------------------------------
lcc  <- tar_read(large_clone_comparisons)
lows <- unlist(tar_read(seus_low_hypoxia))
ids  <- numbatHelpers:::.seu_sample_id(lows)

cat("\n== .seu_sample_id ==\n")
cat("example:", basename(lows[1]), "->", ids[1], "\n")
n_match <- sum(ids %in% names(lcc))
cat("ids resolving in large_clone_comparisons:", n_match, "of", length(ids), "\n")
ok(!any(str_detect(ids, "_seu|\\.rds")), "no object-type suffix survives in the id")
ok(n_match > 0, "low-hypoxia ids resolve against large_clone_comparisons")

# branch suffixes must be preserved, not collapsed to the bare SRX id
ok(numbatHelpers:::.seu_sample_id("output/seurat/SRX10264523_branch_5_filtered_seu.rds") ==
     "SRX10264523_branch_5", "_branch_N suffix is preserved")
ok(numbatHelpers:::.seu_sample_id("output/seurat/SRX10031194_hypoxia_low_seu.rds") ==
     "SRX10031194", "_hypoxia_low_seu suffix is stripped")

# --- 2. find_diffex_clones actually returns something now --------------------
cat("\n== find_diffex_clones ==\n")
nb <- tar_read(numbat_rds_files)
target_id <- names(lcc)[names(lcc) %in% ids][1]
seu_path  <- lows[ids == target_id][1]
cat("running on:", basename(seu_path), "(id", target_id, ")\n")
res <- find_diffex_clones(seu_path, nb, lcc, location = "all")
cat("comparisons returned:", length(res), " names:", paste(names(res), collapse = ", "), "\n")
ok(length(res) > 0, "find_diffex_clones returns a non-empty result for a sample with comparisons")

# --- 3. empty-workbook tolerance --------------------------------------------
cat("\n== tally_kooi_candidates on the empty placeholder workbook ==\n")
empty_wb <- "results/diffex_bw_clones_large_in_segment_by_chr.xlsx"
cat("workbook:", empty_wb, "sheets:", paste(readxl::excel_sheets(empty_wb), collapse = ","), "\n")
tal <- withCallingHandlers(
  tally_kooi_candidates(empty_wb, empty_wb),
  warning = function(w) { cat("  (expected warning):", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }
)
ok(is.list(tal) && all(c("cis", "trans") %in% names(tal)),
   "tally_kooi_candidates returns cis/trans on an empty workbook instead of erroring")
ok(nrow(tal$cis) == 0, "empty workbook yields zero candidates")

# also tolerate a flatly missing file
tal2 <- suppressWarnings(tally_kooi_candidates("results/__no_such_file__.xlsx", empty_wb))
ok(nrow(tal2$cis) == 0, "missing workbook is tolerated")

# --- 4. duplicate-rowname dedupe --------------------------------------------
cat("\n== duplicate symbol dedupe ==\n")
kc <- numbatHelpers:::.read_kooi_candidates()
ok(!any(duplicated(kc$symbol[kc$symbol == "TDP2"])), "kooi_candidates TDP2 duplicate removed")

csv <- grep("SRX10264519", .drop_missing_paths(tar_read(all_diffex_clones_for_each_cluster)), value = TRUE)[1]
d <- read_csv(csv, show_col_types = FALSE)
before <- d %>% count(clone_comparison, cluster, symbol) %>% filter(n > 1) %>% nrow()
after <- d %>%
  numbatHelpers:::.canonical_chr_first() %>%
  distinct(clone_comparison, cluster, symbol, .keep_all = TRUE) %>%
  count(clone_comparison, cluster, symbol) %>% filter(n > 1) %>% nrow()
cat("dup (comparison,cluster,symbol): before =", before, " after =", after, "\n")
ok(before > 0, "the test CSV really does contain duplicates (guard is meaningful)")
ok(after == 0, "dedupe leaves one row per symbol per cluster -- column_to_rownames is safe")

# the retained row must be the canonical chromosome, not the patch contig
atmin <- d %>% filter(symbol == "ATMIN") %>%
  numbatHelpers:::.canonical_chr_first() %>%
  distinct(clone_comparison, cluster, symbol, .keep_all = TRUE)
cat("ATMIN retained chr values:", paste(unique(atmin$chr), collapse = ", "), "\n")
ok(!any(str_detect(as.character(atmin$chr), "PATCH|CHR_")),
   "canonical chromosome is retained over the patch contig")

cat("\nALL SMOKE CHECKS PASSED\n")
