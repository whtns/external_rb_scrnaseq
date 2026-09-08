# Consolidated eligibility + artifact inventory for the RB SCNA clone
# comparisons. Answers: which samples permit a with/without contrast for each
# RB SCNA, of what kind, and what artifacts already exist for them.
#
# Tier A = both clones aneuploid (subclonal SCNA). The tumor background is held
#          fixed, so the contrast isolates the SCNA.
# Tier B = the only SCNA-negative clone is numbat's diploid clone (clonal SCNA).
#          The contrast is tumor-vs-normal and confounds the SCNA with every
#          other difference between tumor and normal cells. These are the
#          samples that need the aggregate diploid object.

suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(stringr)})
setwd("/project2/cobrinik_1090/external_rb_scrnaseq_proj")
out_dir <- "results/diploid_audit"

RB  <- c("1q+", "2p+", "6p+", "16q-")
key <- c("1q+" = "1q", "2p+" = "2p", "6p+" = "6p", "16q-" = "16q")

pr <- read_csv(file.path(out_dir, "clone_pairs_rb_scna.csv"), show_col_types = FALSE) |>
  separate_rows(rb_gained, sep = ";") |>
  filter(rb_gained %in% RB, usable)

pick <- function(d) d |> arrange(n_extra_events, desc(pmin(n_cells_wo, n_cells_w))) |> slice(1)

tierA <- pr |> filter(!wo_is_diploid) |> group_by(sample_id, rb_gained) |> pick() |> ungroup() |> mutate(tier = "A")
tierB <- pr |> filter(wo_is_diploid)  |> group_by(sample_id, rb_gained) |> pick() |> ungroup() |> mutate(tier = "B")
# a sample/SCNA in tier A does not need its tier B entry
tierB <- anti_join(tierB, tierA, by = c("sample_id", "rb_gained"))
elig  <- bind_rows(tierA, tierB) |> mutate(scna = unname(key[rb_gained]))

# current curated collage set
collage <- list(
  "1q"  = c("SRX10264523","SRX10264526","SRX11133594","SRX11133593","SRX11133592","SRX10831287"),
  "2p"  = c("SRX10031193","SRX10264517","SRX10264518","SRX10264519","SRX10264520","SRX10264523","SRX14116946","SRX14116947","SRX22868102"),
  "6p"  = c("SRX10031193","SRX10831281","SRX10831282","SRX11133588","SRX14116944","SRX14116946","SRX14116947","SRX22868105"),
  "16q" = c("SRX11133594","SRX11133593","SRX11133592"))
elig$in_collage_set <- mapply(function(s, k) s %in% collage[[k]], elig$sample_id, elig$scna)

# artifacts on disk
has_file <- function(pat) length(Sys.glob(pat)) > 0
elig$has_summary  <- file.exists(sprintf("results/%s_summary.pdf", elig$sample_id))
elig$n_diffex     <- vapply(elig$sample_id,
  function(s) length(Sys.glob(sprintf("results/%s_*diffex*.csv", s))), integer(1))
elig$has_scna_diffex <- mapply(function(s, r)
  has_file(sprintf("results/%s_*diffex_all_%s.csv", s, r)), elig$sample_id, elig$rb_gained)

elig <- elig |>
  select(sample_id, scna, rb_gained, tier, clone_wo, clone_w, n_cells_wo, n_cells_w,
         n_extra_events, all_gained, in_collage_set, has_summary, n_diffex, has_scna_diffex) |>
  arrange(scna, tier, n_extra_events, desc(pmin(n_cells_wo, n_cells_w)))
write_csv(elig, file.path(out_dir, "rb_scna_comparison_eligibility.csv"))

cat("=== TIER A (tumor-vs-tumor) by SCNA ===\n")
for (k in c("1q","2p","6p","16q")) {
  d <- elig |> filter(scna == k, tier == "A")
  cat("\n--", k, ":", nrow(d), "samples --\n")
  if (nrow(d)) print(as.data.frame(d |> select(sample_id, clone_wo, clone_w, n_cells_wo,
      n_cells_w, n_extra_events, in_collage_set, has_scna_diffex)))
}
cat("\n\n=== collage-set mismatch ===\n")
for (k in c("1q","2p","6p","16q")) {
  a <- elig |> filter(scna == k, tier == "A") |> pull(sample_id)
  cur <- collage[[k]]
  cat("\n", k, ":\n", sep = "")
  cat("  tier A missing from collage set : ", paste(setdiff(a, cur), collapse = ", "), "\n", sep = "")
  cat("  in collage set but NOT tier A   : ", paste(setdiff(cur, a), collapse = ", "), "\n", sep = "")
}
cat("\n=== artifact gaps for tier A ===\n")
g <- elig |> filter(tier == "A") |> group_by(sample_id) |>
  summarise(scnas = paste(rb_gained, collapse=";"), has_summary = all(has_summary),
            n_diffex = max(n_diffex), .groups="drop") |> arrange(n_diffex)
print(as.data.frame(g))
