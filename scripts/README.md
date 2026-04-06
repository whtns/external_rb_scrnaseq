# Scripts Directory

This directory contains utility scripts for the project.

## monitor_targets.sh

A real-time monitoring dashboard for targets pipelines.

### Usage

```bash
# Monitor default store with 10-second refresh
./scripts/monitor_targets.sh

# Monitor specific store
./scripts/monitor_targets.sh _targets

# Custom store and refresh interval (5 seconds)
./scripts/monitor_targets.sh _targets_r431 5
```

### Parameters

1. **store_path** (optional, default: `_targets_r431`)
   - Path to the targets store directory to monitor

2. **refresh_interval** (optional, default: `10`)
   - Seconds between status updates

### Features

- **Progress tracking**: Shows percentage complete and target counts
- **Active targets**: Lists targets currently being executed
- **Error detection**: Alerts if any targets have failed
- **Runtime stats**: Displays recently completed targets with execution times
- **Auto-refresh**: Updates continuously until stopped with Ctrl+C

### Output

The monitor displays:
- Overall progress percentage
- ✓ Completed targets count
- ⊙ Dispatched targets (currently running)
- ⊘ Skipped targets (already up-to-date)
- List of active target names
- Top 5 recently completed targets by runtime
- Error warnings if present

### Environment Variables

- **R_PATH**: Path to Rscript binary (default: `/opt/R/4.3.1/bin/Rscript`)

### Example Output

```
=== TARGETS PIPELINE MONITOR ===
Fri 03 Apr 2026 02:34:11 PM PDT

Store: _targets_r431

Progress: 98.3 %
  ✓ Completed: 22 
  ⊙ Dispatched: 5 
  ⊘ Skipped: 260 
  Total: 287 

Currently active targets:
  • unfiltered_clone_trees_segments_files 
  • unfiltered_clone_trees_segments_files_b5da3d85c932872b 

Recently completed (top 5 by runtime):
  • filtered_seus - 109.3 min 
  • unfiltered_seus - 68.5 min 
  • ideogram_res_s06a - 57.7 min 
```
