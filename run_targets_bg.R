#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
  library(callr)
})

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", timestamp()), ..., "\n", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)

# Any arg starting with "_targets" is the store path; all others are target names.
is_store <- grepl("^_targets", args)
store_path <- if (any(is_store)) args[is_store][[1]] else "_targets_r431"
target_names <- args[!is_store]
if (length(target_names) == 0) target_names <- "figures_and_tables"

log_msg("Starting targets run")
log_msg("Targets: ", paste(target_names, collapse = ", "))
log_msg("Store: ", store_path)
log_msg("Working directory: ", getwd())

runner <- callr::r_bg(
  func = function(target_names, store_path) {
    suppressPackageStartupMessages(library(targets))
    tar_config_set(store = store_path)
    tar_make(names = all_of(target_names), reporter = "timestamp")
  },
  args = list(target_names, store_path),
  stdout = "|",
  stderr = "2>&1",
  supervise = TRUE
)

log_msg("Spawned background R process (pid ", runner$get_pid(), ")")

last_heartbeat <- Sys.time()
heartbeat_seconds <- 30

while (runner$is_alive()) {
  output_lines <- runner$read_output_lines()
  if (length(output_lines)) {
    cat(paste0(output_lines, "\n"), sep = "")
  }

  if (as.numeric(difftime(Sys.time(), last_heartbeat, units = "secs")) >= heartbeat_seconds) {
    log_msg("Still running...")
    last_heartbeat <- Sys.time()
  }

  Sys.sleep(1)
}

# Flush any remaining buffered output.
output_lines <- runner$read_output_lines()
if (length(output_lines)) {
  cat(paste0(output_lines, "\n"), sep = "")
}

exit_status <- runner$get_exit_status()

if (!is.null(exit_status) && exit_status != 0) {
  log_msg("targets run failed with exit status ", exit_status)
  quit(status = exit_status)
}

tar_config_set(store = store_path)
meta <- tar_meta(fields = c(name, error), complete_only = FALSE)
errored <- sum(!is.na(meta$error) & nzchar(meta$error))

if (errored > 0) {
  log_msg("Completed with ", errored, " errored target(s)")
} else {
  log_msg("Completed with 0 errored targets")
}

log_msg("Done")
