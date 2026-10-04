#!/usr/bin/env Rscript
# Joint genotyping + phasing of technical-replicate captures from one tumor.
#
# This is the post-pileup half of numbat's bin/pileup_and_phase.R in
# multi-sample mode (--label <tumor> --samples a,b). The cellsnp-lite pileups
# already exist per replicate under output/numbat/<SRX>/pileup/<SRX>/, so they
# are read in place instead of re-piled. What changes versus the per-replicate
# runs is that genotype() and Eagle see the pooled SNPs once, so both
# replicates' allele tables share one set of phased haplotypes. Rbinding the
# old per-replicate tables would mix independently phased blocks.
#
# Run inside pipeline/containers/numbat-pipeline.sif (numbat 1.5.2, eagle) --
# the same stack that produced output/numbat/<SRX>_allele_counts.tsv.gz.
#
# Usage: Rscript src/phase_merged_replicates.R <label> <SRX_a,SRX_b> <outdir> [ncores]
# Writes: <outdir>/phasing/<label>_chr*.phased.vcf.gz
#         <outdir>/<SRX>_allele_counts.tsv.gz  (barcodes NOT prefixed; see
#         src/build_merged_replicate_inputs.R)

suppressPackageStartupMessages({
  library(numbat)
  library(data.table)
  library(dplyr)
  library(glue)
  library(Matrix)
  library(stringr)
  library(vcfR) # numbat:::genotype() calls write.vcf() unqualified
})

args <- commandArgs(trailingOnly = TRUE)
label <- args[[1]]
samples <- str_split(args[[2]], ",")[[1]]
outdir <- args[[3]]
ncores <- if (length(args) >= 4) as.integer(args[[4]]) else 8

proj <- "/project2/cobrinik_1090/external_rb_scrnaseq_proj"
gmap <- "/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
paneldir <- "/project2/cobrinik_1090/Homo_sapiens/numbat/1000G_hg38"
pileup_dir <- function(s) glue("{proj}/output/numbat/{s}/pileup/{s}")

for (s in samples) {
  stopifnot(file.exists(glue("{pileup_dir(s)}/cellSNP.base.vcf")))
}
dir.create(glue("{outdir}/phasing"), recursive = TRUE, showWarnings = FALSE)

## VCF creation (joint)
cat("Creating VCFs\n")
vcfs <- lapply(samples, function(s) {
  vcf <- vcfR::read.vcfR(glue("{pileup_dir(s)}/cellSNP.base.vcf"), verbose = FALSE)
  if (nrow(vcf@fix) == 0) stop(glue("Pileup VCF for sample {s} has 0 variants"))
  vcf@fix[, 1] <- gsub("chr", "", vcf@fix[, 1])
  vcf
})
numbat:::genotype(label, samples, vcfs, glue("{outdir}/phasing"), chr_prefix = TRUE)

## phasing
cat("Running phasing\n")
cmds <- vapply(1:22, function(chr) {
  paste("eagle",
    glue("--numThreads {ncores}"),
    glue("--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz"),
    glue("--vcfRef {paneldir}/chr{chr}.genotypes.bcf"),
    glue("--geneticMapFile={gmap}"),
    glue("--outPrefix {outdir}/phasing/{label}_chr{chr}.phased"))
}, character(1))
script <- glue("{outdir}/run_phasing.sh")
writeLines(cmds, script)
status <- system2("sh", script, stdout = glue("{outdir}/phasing.log"),
                  stderr = glue("{outdir}/phasing.log"))
if (status != 0) stop("Phasing failed; see phasing.log")

## allele count dataframes, one per replicate, against the shared phasing
cat("Generating allele count dataframes\n")
genetic_map <- fread(gmap) %>%
  setNames(c("CHROM", "POS", "rate", "cM")) %>%
  group_by(CHROM) %>%
  mutate(start = POS, end = c(POS[2:length(POS)], POS[length(POS)])) %>%
  ungroup()

vcf_phased <- lapply(1:22, function(chr) {
  f <- glue("{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz")
  if (!file.exists(f)) stop(glue("Phased VCF not found: {f}"))
  fread(f, skip = "#CHROM") %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, "chr"))
}) %>%
  Reduce(rbind, .) %>%
  mutate(CHROM = factor(CHROM, unique(CHROM)))

for (s in samples) {
  pu <- pileup_dir(s)
  vcf_pu <- fread(glue("{pu}/cellSNP.base.vcf"), skip = "#CHROM") %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, "chr"))
  df <- numbat:::preprocess_allele(
    sample = label,
    vcf_pu = vcf_pu,
    vcf_phased = vcf_phased,
    AD = readMM(glue("{pu}/cellSNP.tag.AD.mtx")),
    DP = readMM(glue("{pu}/cellSNP.tag.DP.mtx")),
    barcodes = fread(glue("{pu}/cellSNP.samples.tsv"), header = FALSE)$V1,
    gtf = gtf_hg38,
    gmap = genetic_map
  ) %>%
    filter(GT %in% c("1|0", "0|1"))
  fwrite(df, glue("{outdir}/{s}_allele_counts.tsv.gz"), sep = "\t")
  cat(s, ": ", nrow(df), " rows, ", n_distinct(df$snp_id), " het SNPs\n", sep = "")
}
cat("All done!\n")
