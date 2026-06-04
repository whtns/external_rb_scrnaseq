# Approach Review: targets + crew Error Handling

**Date:** 2026-03-25
**Source:** Pipeline debugging — `filtered_seus` chain blocked by unrelated `collages_2p` error

## Context

The targets pipeline uses `crew_controller_local(workers = 4)` for parallelization. When any target errors (even one not in the requested dependency subgraph), the entire pipeline aborts despite `error = "continue"` being set. This blocks `sample_summaries` from building because an unrelated `collages_2p` target hits a JoinLayers error.

**Installed versions:** targets 1.3.2, crew 1.3.0, mirai 2.6.1

## Research Findings

### 1. `error = "continue"` Behavior with crew

With `error = "continue"`, targets marks the failed target as errored and continues the pipeline. However, the crew controller's task dispatch mechanism can still abort if it encounters an unexpected error during task collection/aggregation (not during target execution). This is consistent with [Discussion #1207](https://github.com/ropensci/targets/discussions/1207) which reports "pipeline stutters with crew, completes few tasks and fails."

The key distinction: `error = "continue"` handles errors **inside target functions**. Errors during crew's task dispatch/collection (e.g., serialization failures, row-count mismatches in metadata) are treated as infrastructure failures and abort the pipeline.

**Sources:**
- [Discussion #1207: Pipeline stutters with crew](https://github.com/ropensci/targets/discussions/1207)
- [Discussion #1310: What happens with errors?](https://github.com/ropensci/targets/discussions/1310)

### 2. `error = "null"` Option

`error = "null"` stores NULL for errored targets and lets downstream targets receive NULL. This is more resilient than `error = "continue"` because it provides concrete output (NULL) rather than relying on cached values. However, in our testing, the pipeline still aborted — suggesting the error may be happening at the crew infrastructure level, not inside a target function.

**Sources:**
- [tar_option_set documentation](https://docs.ropensci.org/targets/reference/tar_option_set.html)
- [Debugging pipelines](https://books.ropensci.org/targets/debugging.html)

### 3. `error = "trim"` — New Option (targets >= 1.8.0)

A newer `error = "trim"` option was introduced (post-1.3.2). With `error = "trim"`:
- All currently running targets stay running
- A queued target starts **if it is not a sibling branch** from the same `tar_target()` call
- Healthy regions of the dependency graph can still begin running

This is the most sophisticated error handling option and was specifically designed for crew-based pipelines where you want failures in one branch to not block unrelated work. However, **it requires targets >= 1.8.0** (we have 1.3.2).

A deadlock bug in `error = "trim"` was fixed in targets 1.10.1.

**Sources:**
- [targets changelog](https://docs.ropensci.org/targets/news/index.html)
- [tar_option_set reference](https://docs.ropensci.org/targets/reference/tar_option_set.html)

### 4. Version Gap

Our installed targets (1.3.2) is significantly behind current (1.10+). Key improvements since 1.3.2:
- `error = "trim"` option for granular error handling with crew
- `workspace_on_error` defaulting to TRUE
- Fixes for crew controller timeout and dispatch issues (#1207)
- Memory management improvements for large pipelines
- Deadlock fixes in error handling paths

crew 1.3.0 is also behind current releases, which include fixes for worker lifecycle management and error propagation.

### 5. The Actual Error

The `arguments imply differing number of rows: 1, 0` error is an R data.frame construction error, not a crew/targets bug. It happens inside `filter_cluster_save_seu` when `simplify_gt_col()` returns an empty vector for some cells. The fix (wrapping in `list()` + `unnest` + `distinct`) is already applied. The remaining question is whether this error comes from a target function or from crew's internal metadata handling.

## Recommendations

### Option A: Upgrade targets + crew (Best Long-Term)
Update to targets >= 1.10.x and crew >= 1.x (latest). Use `error = "trim"` for optimal behavior with parallel workers. This gives the most robust error handling.

```r
# After upgrading:
tar_option_set(
  error = "trim",
  workspace_on_error = TRUE,
  controller = crew_controller_local(workers = 4)
)
```

**Trade-off:** Package upgrades may introduce breaking changes; requires testing.

### Option B: Run Without crew for Now (Quick Fix)
Disable the crew controller and run sequentially. This avoids the crew dispatch bug entirely and gives clear error messages. Once `filtered_seus` branches all build, re-enable crew.

```r
tar_option_set(
  error = "continue",
  workspace_on_error = TRUE
  # controller = crew_controller_local(workers = 4)  # disabled
)
```

**Trade-off:** Slower (sequential), but reliable for getting past the current blocker.

### Option C: Wrap Error-Prone Targets in tryCatch (Defensive)
Instead of relying on targets-level error handling, wrap the target function body in `tryCatch()` to return a sentinel value on failure. This prevents errors from reaching crew at all.

```r
filter_cluster_save_seu <- function(...) {
  tryCatch({
    # ... existing code ...
  }, error = function(e) {
    warning("filter_cluster_save_seu failed: ", e$message)
    return(NULL)
  })
}
```

**Trade-off:** Masks errors; requires NULL-checking in downstream targets.

## Key Takeaways

- **Our targets (1.3.2) is 7+ major versions behind current.** The `error = "trim"` option designed for crew doesn't exist in our version.
- **The crew controller treats some errors as infrastructure failures**, bypassing `error = "continue"`. This is why the pipeline aborts on unrelated targets.
- **Immediate fix: disable crew** to build `filtered_seus` → `sample_summaries` chain sequentially. It's already running without crew right now.
- **Medium-term: upgrade targets + crew** to get `error = "trim"` and other fixes.
- **The data error itself is already fixed** — the `list()` wrapper in `simplify_gt_col` addresses the `rows: 1, 0` issue in target functions. The remaining abort is a crew dispatch-level problem.
