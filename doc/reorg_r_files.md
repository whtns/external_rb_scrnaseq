Proposed Consolidation: 53 plot_functions_*.R → 9 themed files
Target file	Source plot_functions_*.R	Theme
figure_functions.R	13, 14, 15, 16, 18, 25, 30	plot_fig_* publication/multi-panel figures
clone_functions.R	1, 2, 10, 21, 22, 25, 26, 32, 47, 48, 52	Clone trees, distributions, phylo
scna_functions.R	12, 16, 17, 18, 19, 34, 36, 44, 49, 52	SCNA/CNV visualization & arm analysis
cellcycle_functions.R	8, 20, 43, 50, 52	CC scoring, phase assignment
qc_functions.R	3, 40, 41, 42, 43, 44, 50	QC metrics, filtering, regression effects
marker_functions.R	4, 6, 7, 17, 24, 29, 36, 37, 46, 52, 53	Markers, heatmaps, dotplots
diffex_functions.R (expand existing)	11, 17, 25, 26, 27, 30, 31, 35, 36, 48	Diffex at clone/cluster level
enrichment_functions.R (expand existing)	38, 39, 51	Enrichment analysis & visualization
utility_functions.R	5, 6, 9, 20, 21, 23, 28, 29, 33, 41, 44, 45, 46, 49, 50	Helpers, data wrangling
The existing numbat_functions_*.R, metadata_functions_1.R, integration_functions.R, heatmap_functions.R, seu_metadata_db.R, and send_notification.R are already well-organized and would stay.

This is a large file reorganization (~250+ functions across 53 files). Do you want me to proceed with implementing it? If so, should I:

Do it all at once (move functions, delete old files)?
Start with one theme at a time so you can review each batch?