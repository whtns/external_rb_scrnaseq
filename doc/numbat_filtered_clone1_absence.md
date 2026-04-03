# Why Clone 1 Is Absent from `filtered_clone_tree_files` for Some Samples (e.g. SRR14800541)

## Background

In numbat, clone 1 is the diploid/normal reference clone. It is always **inferred** as a possible state, but it is not guaranteed that any cell will be **assigned** to it as `clone_opt`.

## Investigation Summary

The absence of clone 1 in `filtered_clone_tree_files` for SRR14800541 is **expected and correct**. It is not a bug in `filter_cluster_save_seu` or `plot_clone_tree`.

## Pipeline Flow

1. **`numbat_sridhar_filtered` Snakemake rule** ([pipeline/Snakefile](../pipeline/Snakefile)) runs `scripts/run_numbat_filtered.R` on the pre-filtered Seurat object (`output/seurat/{sample}_filtered_seu.rds`). The count matrix is subsetted to only cells present in that Seurat object before numbat runs — so cells are filtered **before** numbat inference.

2. **`numbat_sridhar_filtered_rds`** then runs `scripts/process_numbat_rds.R` to serialize the numbat output as an RDS.

3. **`filtered_seus_nb_filtered`** calls `filter_cluster_save_seu()` ([R/plot_functions_3.R](../R/plot_functions_3.R)), which reads `clone_opt` and `GT_opt` from the filtered numbat RDS and adds them to the Seurat object. `clone_opt` is **not overwritten** at any later point in that function.

4. **`filtered_clone_tree_files`** calls `save_clone_tree_from_path()` → `plot_clone_tree()` ([R/plot_functions_48.R](../R/plot_functions_48.R)), which filters the mutation graph to `clone %in% unique(seu$clone_opt)`. If no cell has `clone_opt == 1`, clone 1 does not appear in the tree.

## Root Cause

For SRR14800541, the `clone_post_2.tsv` output from the filtered numbat run has columns `p_1` through `p_7`, confirming clone 1 **was inferred**. However, no cell has a maximum posterior probability for clone 1 — every cell's `clone_opt` is 2–7.

This is biologically interpretable: the QC/filtering step that produced `{sample}_filtered_seu.rds` removed all normal (diploid) cells from SRR14800541. With no normal cells in the input, numbat assigns zero cells to clone 1 as their most probable state, even though clone 1 remains a valid latent state in the model.

## Conclusion

- `clone_opt` in the filtered Seurat object is set correctly.
- Clone 1 is absent because the filtered cell population for SRR14800541 is entirely tumor — no diploid cells remain after QC filtering.
- This is expected behavior, not a pipeline bug.
