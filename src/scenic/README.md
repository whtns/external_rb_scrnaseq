
the current workflow for running [SCENIC](https://www.nature.com/articles/s41596-020-0336-2) is based on this github [repo](https://github.com/aertslab/SCENICprotocol) from the aerts lab

It uses nextflow to execute time-intensive processing followed by downstream analysis using a modifiable jupyter notebook. Input is a loom file, possibly converted from a seurat object with cells of interest.
output of the nextflow pipeline is another loom file. 

A helper script, `run_scenic.sh` is included to ease execution of the nextflow pipeline and standardize file paths for input and output loom files. run_scenic.sh also includes scenic-required reference files located at `~/Homo_sapiens/scenic/` that can be swapped out if needed. 

# Quick Start
1. navigate to project directory
1. clone this repository into the `src` subdirectory
`git clone https://github.com/whtns/scenic_src.git`
1. run `run_scenic.sh` with a command-line argument indicating the path to the desired seurat object 
`bash run_scenic.sh unfiltered`
1. open the jupyter notebook `scenic_downstream_analysis.ipynb` in [jupyter lab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) located on the srt server at `cobrinik-srt.la.ad.chla.org:8888` and modify as needed to explore common regulons. modify paths as needed to load `output/scenic/<dataset>-final.loom 

