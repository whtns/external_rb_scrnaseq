#!/opt/R/4.3.1/bin/Rscript

library(argparse)

parser <- ArgumentParser(description='Run VCF genotyping, Eagle phasing, and allele count generation from existing cellsnp-lite pileup')
parser$add_argument('--label',    type = "character", required = TRUE,  help = "Individual label")
parser$add_argument('--samples',  type = "character", required = TRUE,  help = "Sample names, comma delimited")
parser$add_argument('--gmap',     type = "character", required = TRUE,  help = "Path to Eagle2 genetic map")
parser$add_argument('--eagle',    type = "character", required = FALSE, default = 'eagle', help = "Path to Eagle2 binary")
parser$add_argument('--paneldir', type = "character", required = TRUE,  help = "Directory of phasing reference panel (BCF files)")
parser$add_argument('--outdir',   type = "character", required = TRUE,  help = "Output directory (same as used for pileup)")
parser$add_argument('--ncores',   type = "integer",   required = TRUE,  help = "Number of cores")
args <- parser$parse_args()

suppressPackageStartupMessages({
    library(glue)
    library(stringr)
    library(data.table)
    library(dplyr)
    library(vcfR)
    library(Matrix)
    library(numbat)
})

source("scripts/patch_numbat_dplyr.R")

label   <- args$label
samples <- str_split(args$samples, ',')[[1]]
outdir  <- args$outdir
ncores  <- args$ncores
gmap    <- args$gmap
eagle   <- args$eagle
paneldir <- args$paneldir
genome  <- ifelse(str_detect(args$gmap, 'hg19'), 'hg19', 'hg38')
message(paste0('Using genome version: ', genome))

dir.create(glue('{outdir}/phasing'), showWarnings = FALSE, recursive = TRUE)

## VCF creation
cat('Creating VCFs\n')
vcfs <- lapply(samples, function(sample) {
    vcfR::read.vcfR(glue('{outdir}/pileup/{sample}/cellSNP.base.vcf'), verbose = FALSE)
})

# genotype() expects bare chromosome names ("1","2",...); cellsnp-lite outputs
# "chr1","chr2",... — strip the prefix so make_vcf_chr() comparisons work.
vcfs <- lapply(vcfs, function(vcf) {
    vcf@fix[, 1] <- sub("^chr", "", vcf@fix[, 1])
    vcf
})

numbat:::genotype(label, samples, vcfs, glue('{outdir}/phasing'))

## Eagle phasing
eagle_cmd <- function(chr) {
    paste(eagle,
          glue('--numThreads {ncores}'),
          glue('--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz'),
          glue('--vcfRef {paneldir}/chr{chr}.genotypes.bcf'),
          glue('--geneticMapFile={gmap}'),
          glue('--outPrefix {outdir}/phasing/{label}_chr{chr}.phased'),
          sep = ' ')
}

cmds   <- lapply(1:22, function(chr) { eagle_cmd(chr) })
script <- glue('{outdir}/run_phasing.sh')
list(cmds) %>% fwrite(script, sep = '\n')
system(glue('chmod 777 {script}'))

tryCatch({
    system(glue('sh {script}'), intern = TRUE)
}, warning = function(w) {
    stop('Phasing failed')
})

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

gtf <- if (genome == 'hg19') gtf_hg19 else gtf_hg38

genetic_map <- fread(gmap) %>%
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end   = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

for (sample in samples) {

    vcf_phased <- lapply(1:22, function(chr) {
        vcf_file <- glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
        if (file.exists(vcf_file)) {
            data.table::fread(vcf_file) %>%
                rename(CHROM = `#CHROM`) %>%
                mutate(CHROM = str_remove(CHROM, 'chr'))
        }
    }) %>%
        Reduce(rbind, .) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    pu_dir <- glue('{outdir}/pileup/{sample}')
    vcf_pu <- fread(glue('{pu_dir}/cellSNP.base.vcf')) %>% rename(CHROM = `#CHROM`)
    AD     <- readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP     <- readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))
    cell_barcodes <- fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = FALSE) %>% pull(V1)

    df <- numbat:::preprocess_allele(
        sample     = label,
        vcf_pu     = vcf_pu,
        vcf_phased = vcf_phased,
        AD         = AD,
        DP         = DP,
        barcodes   = cell_barcodes,
        gtf        = gtf,
        gmap       = genetic_map
    ) %>%
        filter(GT %in% c('1|0', '0|1'))

    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
}
