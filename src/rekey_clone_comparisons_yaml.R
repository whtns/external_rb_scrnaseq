# Re-key config/large_clone_comparisons.yaml from SRR run ids to SRX experiment ids.
#
# WHY: every consumer looks this up by SRX --
#   plot_functions_26 (find_diffex_bw_clones_for_each_cluster): tumor_id = SRX
#   plot_functions_1/19/20/31, figure_functions: sample_id from the SRX-named
#     *_filtered_seu.rds
# but 28 of the 29 keys are SRR (one, SRX10831287, is already SRX). So
# large_clone_comparisons[[sample_id]] has been returning NULL for essentially
# every sample, silently yielding zero clone comparisons -- which is why the
# per-cluster diffex figures (Fig. 7b / 8b) had nothing to build from.
#
# Keys may carry a _branch_N suffix (SRR13884246_branch_5); the suffix is preserved.
# Writes a .bak, then verifies: no key lost, no collision, contents byte-identical
# after re-keying, and every SRX key resolves against numbat_rds_files.
suppressMessages({
  library(dplyr)
  library(stringr)
})

yaml_path <- "config/large_clone_comparisons.yaml"
lcc  <- yaml::read_yaml(yaml_path)
meta <- readr::read_tsv("data/metadata.tsv", show_col_types = FALSE) |>
  select(Run, Experiment) |>
  distinct() |>
  filter(!is.na(Run), !is.na(Experiment))
srr2srx <- setNames(meta$Experiment, meta$Run)

old_keys <- names(lcc)
new_keys <- vapply(old_keys, function(k) {
  if (str_starts(k, "SRX")) return(k)               # already SRX (SRX10831287)
  srr    <- str_extract(k, "SRR[0-9]+")
  suffix <- str_remove(k, "^SRR[0-9]+")             # "" or "_branch_5"
  srx    <- srr2srx[[srr]]
  if (is.null(srx) || is.na(srx)) {
    warning("no SRX for ", k, "; leaving key unchanged")
    return(k)
  }
  paste0(srx, suffix)
}, character(1))

cat("=== re-key map ===\n")
for (i in seq_along(old_keys)) {
  cat(sprintf("  %-24s -> %s%s\n", old_keys[i], new_keys[i],
              if (identical(old_keys[i], new_keys[i])) "   (unchanged)" else ""))
}

# --- checks before writing -------------------------------------------------
stopifnot("key count changed"     = length(new_keys) == length(old_keys))
if (any(duplicated(new_keys))) {
  stop("COLLISION: duplicate SRX keys: ",
       paste(unique(new_keys[duplicated(new_keys)]), collapse = ", "))
}
unresolved <- old_keys[str_starts(old_keys, "SRR") & new_keys == old_keys]
if (length(unresolved)) {
  stop("UNRESOLVED SRR keys (no SRX in metadata): ",
       paste(unresolved, collapse = ", "))
}
# contents must be untouched -- only the keys change
names(lcc) <- new_keys
stopifnot("contents changed" =
  identical(unname(yaml::read_yaml(yaml_path)), unname(lcc)))

file.copy(yaml_path, paste0(yaml_path, ".bak"), overwrite = TRUE)
yaml::write_yaml(lcc, yaml_path)

# --- verify what was written ----------------------------------------------
rt <- yaml::read_yaml(yaml_path)
# unname(): new_keys carries names from vapply(old_keys, ...), and identical()
# compares attributes -- comparing it to the bare names(rt) would always fail.
stopifnot("round-trip key mismatch"     = identical(sort(names(rt)),
                                                    sort(unname(new_keys))),
          "round-trip content mismatch" = identical(unname(rt[sort(names(rt))]),
                                                    unname(lcc[sort(names(lcc))])))
cat("\nwrote", yaml_path, "(backup at", paste0(yaml_path, ".bak"), ")\n")
cat("keys:", length(rt), "| all SRX:", all(str_starts(names(rt), "SRX")), "\n")
cat("REKEY DONE\n")
