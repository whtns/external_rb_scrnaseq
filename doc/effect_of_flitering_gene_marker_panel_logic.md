Panel information flow: SRR14800534 effect-of-filtering
Data lineage

unfiltered_seu
    │  all cells, clusters 0–6 computed by FindClusters(gene_snn_res.0.2)
    │  cluster 5 = low-quality / non-tumor island (marked remove=1 in cluster_dictionary)
    │
    ▼ filter_cluster_save_seu()  — removes cluster 5 cells (+ qc_fail, clone_opt_na, etc.)
filtered_seu
    │  clusters 0,1,2,3,4,6  (cluster 5 absent; labels inherited, not re-assigned)
    │
    ▼ load_and_save_hypoxia_score()  — adds hypoxia_score to metadata, no cell removal
hypoxia_seu
    │
    ▼ subset_seu_by_expression(hypoxia_score <= threshold, run_hypoxia_clustering=TRUE)
seus_low_hypoxia
       cells with low hypoxia_score, RE-CLUSTERED at gene_snn_res.0.2
       clusters 0–6  (7 clusters; labels are new assignments, unrelated to unfiltered cluster 5)
Why the cluster counts differ
Panel	Source	Cluster labels	N clusters
unfiltered	unfiltered_seu	original FindClusters	7 (0–6)
filtered	filtered_seus	inherited from unfiltered, cluster 5 dropped	6 (0,1,2,3,4,6)
low_hypoxia	seus_low_hypoxia	re-computed on the hypoxia subset	7 (0–6)
Key point: the cluster 5 visible in the low_hypoxia panel is not the same population as cluster 5 in the unfiltered panel. It is a new cluster produced by re-clustering a smaller, hypoxia-filtered cell set. The shared label is coincidental.