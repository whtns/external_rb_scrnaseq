Now I can see the full loop structure. The answer is a combination of a few things:

1. The loop grows, not shrinks, between iterations

The subtrees variable expands each iteration as the phylogeny resolves more clones. The log shows this directly: iteration 1 runs HMMs on 9 cell groups, iteration 2 runs on 15 cell groups. More groups = more bulk data structures allocated on top of everything else.

2. Critical objects are never freed within the loop

The large inputs — count_mat, lambdas_ref, df_allele, sc_refs — live in the function's enclosing environment and are passed into every call inside the loop. They're never rm()-d or scoped to iteration. So iteration 2 starts with the full baseline footprint already committed.

3. Iteration-scoped objects are reassigned but not promptly collected

Objects like bulk_subtrees, bulk_clones, exp_post, dist_mat, P are reassigned at the top of the next iteration, which makes the old values eligible for GC — but R's garbage collector is lazy: it only runs when the allocator decides it needs memory. There's no explicit gc() call between iterations, so the old iteration 1 objects can be sitting live in memory at the same time new iteration 2 objects are being allocated.

4. The per-cell evaluation step is inherently bigger in iteration 2

get_exp_post (the "Evaluating CNV per cell" step that killed the process) scores every cell against segs_consensus. After iteration 1, segs_consensus has more refined segments than the initial estimate. The memory cost of that step scales as O(cells × segments), so iteration 2 is doing materially more work than iteration 1 did at the same step.

The net effect: iteration 2 inherits the full baseline memory of the function, then allocates new structures proportional to a more complex model, while the GC hasn't had a chance to reclaim iteration 1's temporaries.