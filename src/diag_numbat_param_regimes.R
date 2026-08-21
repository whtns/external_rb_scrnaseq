# What parameters did each production numbat run actually use, and what did each
# regime recover?
#
# Every run writes its full parameter set to output/numbat_sridhar/<sample>/log.txt
# ("Running under parameters:"). Those values are authoritative -- they are what
# numbat itself echoed -- and they disagree with the current pipeline/config.yaml,
# which has drifted (config now says alpha 1e-4; 31 samples ran at 1e-3).
#
# The cohort turns out to span several regimes, including 17 samples at numbat's
# default t = 1e-5 against 54 at t = 1e-2. That makes the t comparison partly
# answerable from data already on disk -- observationally, and confounded by
# study, but informative about whether a t sweep is worth the compute.
#
# Depends on results/rb_scna_canonical_by_sample.csv from
# src/diag_rb_scna_recovery.R.
#
# Output: results/numbat_param_regimes.csv + printed summary.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tidyr)
})

NUMBAT_DIR <- "output/numbat_sridhar"
CANON      <- "results/rb_scna_canonical_by_sample.csv"

# numbat writes "key = value" lines between "Running under parameters:" and
# "Input metrics:".
parse_log <- function(f) {
  lines <- readLines(f, warn = FALSE, n = 60)
  start <- grep("Running under parameters:", lines, fixed = TRUE)
  stop  <- grep("Input metrics:", lines, fixed = TRUE)
  if (length(start) == 0) return(NULL)
  stop <- if (length(stop) == 0) length(lines) else stop[1]
  kv <- lines[(start[1] + 1):(stop - 1)]
  kv <- kv[str_detect(kv, "^\\s*[A-Za-z_]+ = ")]
  tibble(
    sample_id = basename(dirname(f)),
    key   = str_trim(str_match(kv, "^\\s*([A-Za-z_]+) = ")[, 2]),
    value = str_trim(str_match(kv, " = (.*)$")[, 2])
  )
}

logs <- list.files(NUMBAT_DIR, pattern = "^log\\.txt$", recursive = TRUE, full.names = TRUE)
message("parsing ", length(logs), " run logs")

params <- map_dfr(logs, parse_log) %>%
  filter(key %in% c("t", "alpha", "gamma", "min_cells", "init_k", "max_iter",
                    "min_LLR", "max_entropy", "tau", "multi_allelic",
                    "check_convergence", "min_overlap", "skip_nj", "use_loh",
                    "segs_loh", "diploid_chroms", "ncores")) %>%
  pivot_wider(names_from = key, values_from = value)

canon <- read_csv(CANON, show_col_types = FALSE)

recovery <- canon %>%
  group_by(sample_id) %>%
  summarise(
    n_final = sum(called_final),
    n_union = n(),
    .groups = "drop"
  )

tbl <- params %>%
  left_join(recovery, by = "sample_id") %>%
  mutate(across(c(n_final, n_union), ~ replace_na(.x, 0L))) %>%
  arrange(t, alpha, max_iter, sample_id)

dir.create("results", showWarnings = FALSE)
write_csv(tbl, "results/numbat_param_regimes.csv")

cat("\n=== parameter regimes actually run (from log.txt) ===\n")
tbl %>%
  count(t, alpha, gamma, min_cells, max_iter, min_LLR, max_entropy, multi_allelic,
        name = "n_samples") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== canonical RB events recovered, by regime ===\n")
cat("(mean events per sample; max possible = 5)\n\n")
tbl %>%
  group_by(t, alpha, max_iter) %>%
  summarise(
    n_samples     = n(),
    mean_final    = round(mean(n_final), 2),
    mean_union    = round(mean(n_union), 2),
    zero_event    = sum(n_final == 0),
    .groups = "drop"
  ) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== t = 1e-5 (numbat default) vs t = 1e-2 (this pipeline) ===\n")
cat("Confounded by study/sample -- the t = 1e-5 group is not a random subset --\n")
cat("so read this as a prior on whether a controlled sweep is worth running.\n\n")
tbl %>%
  mutate(t_group = ifelse(as.numeric(t) <= 1e-4, "t = 1e-5 (default)", "t = 1e-2 (pipeline)")) %>%
  group_by(t_group) %>%
  summarise(
    n_samples  = n(),
    mean_final = round(mean(n_final), 2),
    mean_union = round(mean(n_union), 2),
    zero_event = sum(n_final == 0),
    .groups = "drop"
  ) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nwrote results/numbat_param_regimes.csv\n")
