# Dockerfile for Snakemake numbat-only workflow
#
# Rule chain (rule all):
#   pileup → phase_and_count → seurat → numbat_sridhar → numbat_sridhar_rds
#
# Build:
#   docker build -t numbat-pipeline .
#
# Run (numbat needs up to ~50 GB RAM; allocate headroom above numbat_mem_mb):
#   docker run --rm \
#     --memory=70g --memory-swap=70g \
#     --ulimit nofile=65536:65536 \
#     -v /path/to/project:/proj \
#     -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
#     -w /proj/pipeline \
#     numbat-pipeline \
#     snakemake --snakefile Snakefile --cores 6
#
# NOTE: the numbat rule enforces its own virtual-memory ceiling via `ulimit -v`
# using the value of numbat_mem_mb in config.yaml (currently 65536 MB = 64 GB).
# --memory must be >= that value or the container OOM-killer fires first.

FROM pkharchenkolab/numbat-rbase:latest

LABEL description="Numbat scRNA-seq CNV pipeline (inference + RDS only)"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC

# ---------------------------------------------------------------------------
# System dependencies (single layer)
#   libglpk-dev libgmp-dev    — igraph (Seurat / numbat)
#   libhdf5-dev               — hdf5r (Seurat I/O)
#   libharfbuzz-dev           — textshaping → ragg → ggplot2
#   libfribidi-dev            — textshaping
#   libfontconfig1-dev        — systemfonts → ggplot2
#   libtbb-dev                — fgsea C++ backend
#   python3-pip python3.13-dev — Snakemake (datrie needs Python.h)
# ---------------------------------------------------------------------------
RUN apt-get update --allow-releaseinfo-change && apt-get install -y --no-install-recommends \
    libglpk-dev \
    libgmp-dev \
    libhdf5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libtbb-dev \
    python3-pip \
    python3.13-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# R packages — CRAN
# Base image (Debian trixie) ships Matrix 1.6-1, but SeuratObject >= 5 requires
# Matrix >= 1.6.4.  CRAN's current Matrix (1.7+) needs R >= 4.4.0 and is
# therefore invisible to available.packages() on R 4.3.1.  Install 1.6-4 from
# the CRAN archive first so the Seurat dependency chain resolves.
# Seurat uses restricted deps to skip Suggests (BPCells, hdf5r, etc.) that
# are unused in this pipeline's Read10X → saveRDS workflow.
# ---------------------------------------------------------------------------
RUN Rscript -e "\
    options(repos = c(CRAN='https://cloud.r-project.org')); \
    install.packages( \
        'https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz', \
        repos=NULL, type='source'); \
    stopifnot(packageVersion('Matrix') >= '1.6.4')"

RUN Rscript -e "\
    options(repos = c(CRAN='https://cloud.r-project.org')); \
    install.packages(c('Seurat', 'SeuratObject'), \
        dependencies=c('Depends','Imports','LinkingTo'), Ncpus=4); \
    install.packages(c('conflicted', 'argparse', 'vcfR', 'tidyverse'), \
        dependencies=TRUE, Ncpus=4); \
    stopifnot(requireNamespace('Seurat',    quietly=TRUE)); \
    stopifnot(requireNamespace('tidyverse', quietly=TRUE))"

# ---------------------------------------------------------------------------
# R packages — Bioconductor
# scuttle: used by run_numbat.R via scuttle::downsampleMatrix()
# fgsea: NOT needed — not in numbat Imports and unused by pipeline scripts
# ---------------------------------------------------------------------------
RUN Rscript -e "\
    if (!requireNamespace('scuttle', quietly=TRUE)) BiocManager::install('scuttle', ask=FALSE, update=FALSE, Ncpus=4); \
    stopifnot(requireNamespace('scuttle', quietly=TRUE))"

# ---------------------------------------------------------------------------
# R packages — GitHub
# presto: fast Wilcoxon marker test (used in process_seurat.R)
# numbat: whtns fork — includes per-iteration rm()/gc() memory patches,
#   seq_len(n()) dplyr fix, and make_vcf_chr CHROM factor comparison fix
#   applied directly in source.
# ---------------------------------------------------------------------------
RUN Rscript -e "\
    remotes::install_github('immunogenomics/presto', upgrade='never'); \
    remotes::install_github('whtns/numbat',          upgrade='never', Ncpus=4); \
    stopifnot(requireNamespace('presto', quietly=TRUE)); \
    stopifnot(requireNamespace('numbat', quietly=TRUE))"

# ---------------------------------------------------------------------------
# Snakemake + pandas via pip
# pulp renamed list_solvers → listSolvers in 2.7; patch snakemake 7.x to match
# ---------------------------------------------------------------------------
RUN python3.13 -m pip install --no-cache-dir --break-system-packages \
    "snakemake==7.32.4" "pandas>=1.3"

RUN python3.13 -c "\
import pathlib; \
p = pathlib.Path('/usr/local/lib/python3.13/dist-packages/snakemake/__init__.py'); \
txt = p.read_text(); \
patched = txt.replace('pulp.list_solvers(', 'pulp.listSolvers('); \
p.write_text(patched); \
print('Patched', txt.count('pulp.list_solvers('), 'occurrence(s)')"

# ---------------------------------------------------------------------------
# Verify every R package required by the pipeline scripts is loadable.
# install.packages() and BiocManager::install() exit 0 on silent failure, so
# this layer is the definitive gate before the image is usable.
# ---------------------------------------------------------------------------
RUN Rscript -e "\
    pkgs <- c('Seurat', 'SeuratObject', 'tidyverse', 'fs', 'presto', \
              'numbat', 'argparse', 'vcfR', 'scuttle'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing) > 0) stop('Missing R packages: ', paste(missing, collapse=', '))"

# ---------------------------------------------------------------------------
# Copy pipeline into image
# ---------------------------------------------------------------------------
WORKDIR /pipeline
COPY . /pipeline/

# ---------------------------------------------------------------------------
# Smoke-test: confirm R, Rscript, and snakemake are reachable
# ---------------------------------------------------------------------------
RUN R --version | head -1 \
    && Rscript --version \
    && snakemake --version

CMD ["snakemake", "--snakefile", "Snakefile", "--cores", "1", "--help"]
