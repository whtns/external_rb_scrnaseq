#!/bin/bash
# Generalized targets pipeline monitoring script
# Usage: ./monitor_targets.sh [store_path] [refresh_interval]
#   store_path: Path to targets store (default: _targets_r431)
#   refresh_interval: Seconds between refreshes (default: 10)

STORE_PATH="${1:-_targets_r431}"
REFRESH_INTERVAL="${2:-10}"
R_PATH="${R_PATH:-/opt/R/4.3.1/bin/Rscript}"

echo "Monitoring targets store: $STORE_PATH"
echo "Refresh interval: ${REFRESH_INTERVAL}s"
echo "Press Ctrl+C to stop monitoring"
echo ""

while true; do
  clear
  echo "=== TARGETS PIPELINE MONITOR ==="
  date
  echo ""
  echo "Store: $STORE_PATH"
  echo ""
  
  $R_PATH -e "
    library(targets)
    tar_config_set(store = '$STORE_PATH')
    
    progress <- tar_progress()
    status_table <- table(progress\$progress)
    
    total <- nrow(progress)
    completed <- sum(progress\$progress == 'completed')
    dispatched <- sum(progress\$progress == 'dispatched')
    skipped <- sum(progress\$progress == 'skipped')
    running <- sum(progress\$progress == 'running')
    
    pct_done <- round((completed + skipped) / total * 100, 1)
    
    cat('Progress:', pct_done, '%\n')
    cat('  ✓ Completed:', completed, '\n')
    cat('  ⊙ Dispatched:', dispatched, '\n')
    cat('  ⊘ Skipped:', skipped, '\n')
    cat('  Total:', total, '\n\n')
    
    if (dispatched > 0 || running > 0) {
      cat('Currently active targets:\n')
      active <- progress[progress\$progress %in% c('dispatched', 'running'), ]
      for(i in 1:nrow(active)) {
        cat('  •', active\$name[i], '\n')
      }
      cat('\n')
    } else if (completed + skipped == total) {
      cat('✓ Pipeline complete!\n\n')
    } else {
      cat('No active targets (pipeline may be idle)\n\n')
    }
    
    # Check for errors
    meta <- tar_meta(fields = c(name, error), complete_only = FALSE)
    errors <- meta[!is.na(meta\$error) & nzchar(meta\$error), ]
    if (nrow(errors) > 0) {
      cat('⚠️  Errors detected:', nrow(errors), '\n')
      cat('   Run tar_meta() to see details\n\n')
    }
    
    # Show recently completed targets
    if (completed > 0) {
      meta_all <- tar_meta(fields = c(name, seconds))
      recent <- meta_all[!is.na(meta_all\$seconds), ]
      if (nrow(recent) > 0) {
        recent <- recent[order(recent\$seconds, decreasing = TRUE), ]
        cat('Recently completed (top 5 by runtime):\n')
        for(i in 1:min(5, nrow(recent))) {
          runtime <- round(recent\$seconds[i], 1)
          if (runtime >= 60) {
            runtime_str <- sprintf('%.1f min', runtime / 60)
          } else {
            runtime_str <- sprintf('%.1f sec', runtime)
          }
          cat('  •', recent\$name[i], '-', runtime_str, '\n')
        }
      }
    }
  " 2>&1
  
  echo ""
  echo "Press Ctrl+C to stop monitoring"
  echo "Refreshing in ${REFRESH_INTERVAL} seconds..."
  sleep $REFRESH_INTERVAL
done
