#!/usr/bin Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

# Rprof(rprof_out)

library(numbat)
library(Seurat)
library(readr)
library(magrittr)
library(fs)
conflicted::conflict_prefer("rowSums", "Matrix")

numbat_dir = fs::path_dir(done_file)

## ------------------------------------------------------------------------------------ ##
## Which consensus round to build the object from.
##
## Numbat$new() defaults to i = 2. That is not the final round and never was a
## deliberate choice here -- it is simply the package default, and this script
## previously called Numbat$new(out_dir) with no i, so every RDS in the cohort
## was built from round 2 regardless of how many rounds the sample ran.
##
## The consensus iteration is not monotone (canonical RB events per round run
## 106/98/92/86 across the 33 clean 4-round samples), so no fixed i is right for
## every sample. src/select_numbat_round.R scores each round per sample and
## records a favorable one; results/numbat_selected_round.csv is that manifest.
##
## Behaviour: if the manifest exists and lists this sample, use its round.
## Otherwise fall back to the previous default so unlisted samples are unchanged.
## ------------------------------------------------------------------------------------ ##

DEFAULT_ITER <- 2L

resolve_selected_iter <- function(out_dir) {
  sid <- basename(out_dir)

  # out_dir is <root>/output/<numbat_dirname>/<sample>; the manifest lives at
  # <root>/results/. Derived from out_dir so this works regardless of the
  # working directory snakemake happens to invoke us from.
  root <- fs::path_dir(fs::path_dir(fs::path_dir(out_dir)))
  manifest <- if (exists("selection_manifest")) {
    selection_manifest
  } else {
    fs::path(root, "results", "numbat_selected_round.csv")
  }

  if (!fs::file_exists(manifest)) {
    message(sprintf("[round] no manifest at %s; using default i = %d",
                    manifest, DEFAULT_ITER))
    return(DEFAULT_ITER)
  }

  m <- tryCatch(data.table::fread(manifest), error = function(e) NULL)
  if (is.null(m) || !all(c("sample_id", "selected_round") %in% names(m))) {
    warning(sprintf("[round] manifest %s unreadable; using default i = %d",
                    manifest, DEFAULT_ITER))
    return(DEFAULT_ITER)
  }

  hit <- m[m$sample_id == sid, ]
  if (nrow(hit) != 1L) {
    warning(sprintf("[round] %s not in manifest; using default i = %d",
                    sid, DEFAULT_ITER))
    return(DEFAULT_ITER)
  }

  k <- as.integer(hit$selected_round[1])

  # Only honour the manifest if the round's artefacts are actually on disk.
  need <- c(sprintf("segs_consensus_%d.tsv", k), sprintf("clone_post_%d.tsv", k),
            sprintf("joint_post_%d.tsv", k), sprintf("exp_post_%d.tsv", k),
            sprintf("allele_post_%d.tsv", k), sprintf("geno_%d.tsv", k),
            sprintf("treeML_%d.rds", k), sprintf("mut_graph_%d.rds", k),
            sprintf("tree_final_%d.rds", k))
  missing <- need[!fs::file_exists(fs::path(out_dir, need))]
  if (length(missing) > 0) {
    warning(sprintf("[round] %s selected round %d but missing %s; using default i = %d",
                    sid, k, paste(missing, collapse = ", "), DEFAULT_ITER))
    return(DEFAULT_ITER)
  }

  message(sprintf("[round] %s using selected round i = %d", sid, k))
  k
}

target_i <- resolve_selected_iter(numbat_dir)

# Prefer the highest existing iteration on disk; fall back to max_iter from log.
detect_target_iter <- function(out_dir) {
  iter_files <- fs::dir_ls(out_dir, regexp = "_(\\d+)\\.(tsv|tsv\\.gz)$", type = "file")
  if (length(iter_files) > 0) {
    iters <- as.integer(sub(".*_([0-9]+)\\.(tsv|tsv\\.gz)$", "\\1", basename(iter_files)))
    iters <- iters[!is.na(iters)]
    if (length(iters) > 0) {
      return(max(iters))
    }
  }

  log_path <- fs::path(out_dir, "log.txt")
  if (fs::file_exists(log_path)) {
    lines <- readLines(log_path, warn = FALSE)
    hit <- grep("max_iter\\s*=", lines, value = TRUE)
    if (length(hit) > 0) {
      val <- sub(".*max_iter\\s*=\\s*([0-9]+).*", "\\1", hit[1])
      if (grepl("^[0-9]+$", val)) {
        return(as.integer(val))
      }
    }
  }

  1L
}

# Create missing files expected at target_iter by copying from the highest available
# iteration-specific file for each prefix.
ensure_iteration_aliases <- function(out_dir, target_iter) {
  prefixes <- c("joint_post", "exp_post", "allele_post", "segs_consensus", "geno")
  all_files <- fs::dir_ls(out_dir, type = "file")
  all_basenames <- basename(all_files)

  for (prefix in prefixes) {
    keep <- grepl(paste0("^", prefix, "_(\\d+)\\.tsv$"), all_basenames)
    matches <- all_files[keep]
    if (length(matches) == 0) {
      next
    }

    iters <- as.integer(sub(paste0("^", prefix, "_([0-9]+)\\.tsv$"), "\\1", basename(matches)))
    valid <- which(!is.na(iters))
    if (length(valid) == 0) {
      next
    }

    matches <- matches[valid]
    iters <- iters[valid]
    source_file <- matches[which.max(iters)]
    target_file <- fs::path(out_dir, paste0(prefix, "_", target_iter, ".tsv"))

    if (!fs::file_exists(target_file)) {
      file.copy(source_file, target_file, overwrite = FALSE)
    }
  }

  # numbat sometimes expects this non-iterated final table.
  final_bulk <- fs::path(out_dir, "bulk_clones_final.tsv.gz")
  if (!fs::file_exists(final_bulk)) {
    bulk_keep <- grepl("^bulk_clones_(\\d+)\\.tsv\\.gz$", all_basenames)
    bulk_matches <- all_files[bulk_keep]
    if (length(bulk_matches) > 0) {
      bulk_iters <- as.integer(sub("^bulk_clones_([0-9]+)\\.tsv\\.gz$", "\\1", basename(bulk_matches)))
      valid <- which(!is.na(bulk_iters))
      if (length(valid) > 0) {
        bulk_matches <- bulk_matches[valid]
        bulk_iters <- bulk_iters[valid]
        bulk_source <- bulk_matches[which.max(bulk_iters)]
        file.copy(bulk_source, final_bulk, overwrite = FALSE)
      }
    }
  }
}

build_fallback_nb <- function(out_dir, target_iter) {
  read_if_exists <- function(path, reader) {
    if (fs::file_exists(path)) {
      return(reader(path))
    }
    NULL
  }

  fallback <- list(
    out_dir = out_dir,
    recovered = TRUE,
    recovered_iter = target_iter,
    clone_post = read_if_exists(
      fs::path(out_dir, paste0("joint_post_", target_iter, ".tsv")),
      data.table::fread
    ),
    joint_post = read_if_exists(
      fs::path(out_dir, paste0("joint_post_", target_iter, ".tsv")),
      data.table::fread
    ),
    exp_post = read_if_exists(
      fs::path(out_dir, paste0("exp_post_", target_iter, ".tsv")),
      data.table::fread
    ),
    allele_post = read_if_exists(
      fs::path(out_dir, paste0("allele_post_", target_iter, ".tsv")),
      data.table::fread
    ),
    segs_consensus = read_if_exists(
      fs::path(out_dir, paste0("segs_consensus_", target_iter, ".tsv")),
      data.table::fread
    ),
    bulk_clones = read_if_exists(
      fs::path(out_dir, "bulk_clones_final.tsv.gz"),
      data.table::fread
    )
  )

  class(fallback) <- c("numbat_recovered", "list")
  fallback
}

nb <- tryCatch(
  Numbat$new(out_dir = numbat_dir, i = target_i),
  error = function(e) {
    target_iter <- detect_target_iter(numbat_dir)
    ensure_iteration_aliases(numbat_dir, target_iter)
    tryCatch(
      Numbat$new(out_dir = numbat_dir, i = target_iter),
      error = function(e2) {
        warning(sprintf(
          "Numbat$new failed after recovery for %s; saving fallback object from iteration %d.",
          numbat_dir,
          target_iter
        ))
        build_fallback_nb(numbat_dir, target_iter)
      }
    )
  }
)

## Numbat$new always reads bulk_clones_final.tsv.gz, which comes from the LAST
## round, so on any sample where the selected round is not the last the object
## would mix a round-k segmentation with round-final pseudobulk. Repoint it at
## the selected round's own bulk when that file exists, so the object is
## internally consistent.
bulk_k <- fs::path(numbat_dir, sprintf("bulk_clones_%d.tsv.gz", target_i))
if (inherits(nb, "R6") && fs::file_exists(bulk_k)) {
  ok <- tryCatch({
    nb$bulk_clones <- numbat:::relevel_chrom(data.table::fread(bulk_k))
    TRUE
  }, error = function(e) FALSE)
  message(sprintf("[round] bulk_clones repointed to round %d: %s", target_i, ok))
}

# Record provenance on the object where the class allows it.
try(nb$selected_round <- target_i, silent = TRUE)

nb_path = paste0(fs::path(numbat_dir), "_numbat.rds")

saveRDS(nb, nb_path)