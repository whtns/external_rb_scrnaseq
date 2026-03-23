#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript export_targets_rulegraph.R [output_html] [store_path] [targets_only] [arrangement]",
      "",
      "Arguments:",
      "  output_html  Output HTML path (default: doc/targets_rulegraph.html)",
      "  store_path   targets store path (default: _targets_r431)",
      "  targets_only TRUE/FALSE to hide helper nodes (default: TRUE)",
      "  arrangement  one of: hierarchical-lr, hierarchical-ud, force (default: hierarchical-lr)",
      sep = "\n"
    ),
    "\n"
  )
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  x <- tolower(x)
  if (x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop("Invalid boolean value for targets_only: ", x)
}

apply_arrangement <- function(graph, arrangement) {
  arrangement <- tolower(arrangement)

  if (arrangement == "hierarchical-lr") {
    return(
      graph |>
        visNetwork::visHierarchicalLayout(
          enabled = TRUE,
          direction = "LR",
          sortMethod = "directed",
          levelSeparation = 260,
          nodeSpacing = 220,
          treeSpacing = 320
        ) |>
        visNetwork::visPhysics(enabled = FALSE)
    )
  }

  if (arrangement == "hierarchical-ud") {
    return(
      graph |>
        visNetwork::visHierarchicalLayout(
          enabled = TRUE,
          direction = "UD",
          sortMethod = "directed",
          levelSeparation = 210,
          nodeSpacing = 180,
          treeSpacing = 240
        ) |>
        visNetwork::visPhysics(enabled = FALSE)
    )
  }

  if (arrangement == "force") {
    return(
      graph |>
        visNetwork::visPhysics(
          enabled = TRUE,
          solver = "forceAtlas2Based",
          stabilization = list(enabled = TRUE, iterations = 250)
        )
    )
  }

  stop(
    "Invalid arrangement: ", arrangement,
    ". Use one of hierarchical-lr, hierarchical-ud, force"
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[[1]] %in% c("-h", "--help")) {
  usage()
  quit(status = 0)
}

output_html <- if (length(args) >= 1) args[[1]] else "doc/targets_rulegraph.html"
store_path <- if (length(args) >= 2) args[[2]] else "_targets_r431"
targets_only <- parse_bool(if (length(args) >= 3) args[[3]] else NULL, default = TRUE)
arrangement <- if (length(args) >= 4) args[[4]] else "hierarchical-lr"

if (!dir.exists(store_path)) {
  stop("Store path does not exist: ", store_path)
}

if (!requireNamespace("visNetwork", quietly = TRUE)) {
  stop(
    "Package 'visNetwork' is required to save the graph. ",
    "Install with: install.packages('visNetwork')"
  )
}

out_dir <- dirname(output_html)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

tar_config_set(store = store_path)
graph <- tar_visnetwork(targets_only = targets_only)
graph <- apply_arrangement(graph, arrangement)
graph <- graph |>
  visNetwork::visNodes(font = list(size = 22)) |>
  visNetwork::visEdges(arrows = "to", smooth = FALSE)
visNetwork::visSave(graph, file = output_html)

cat("Saved targets rulegraph to:", output_html, "\n")
cat("Using store:", store_path, "\n")
cat("targets_only:", targets_only, "\n")
cat("arrangement:", arrangement, "\n")