library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## Get the TwoBit with the genomic sequence matching the Ensembl version
## using the AnnotationHub package.
dna <- ensembldb:::getGenomeTwoBitFile(edb)

## Get start/end coordinates of all genes.
genes <- genes(edb)
## Subset to all genes that are encoded on chromosomes for which
## we do have DNA sequence available.
genes <- genes[seqnames(genes) %in% seqnames(seqinfo(dna))]

thrb <- genes[genes$symbol == "THRB"]

## Get start/end coordinates of all transcripts.
transcripts <- transcripts(edb)
## Subset to all transcripts that are encoded on chromosomes for which
## we do have DNA sequence available.
transcripts <- transcripts[seqnames(transcripts) %in% seqnames(seqinfo(dna))]

thrb_tx <-
  transcripts %>%
  as_tibble() %>%
  dplyr::filter(gene_id == thrb$gene_id) %>%
  dplyr::filter(tx_id == "ENST00000396671") %>%
  identity()

## get all exons of all transcripts encoded on chromosome Y
thrbTx <- exonsBy(edb, filter = TxIdFilter("ENST00000396671"))

## Retrieve the sequences for these transcripts from the TwoBitile.
library(GenomicFeatures)
thrbTxSeqs <- extractTranscriptSeqs(dna, thrbTx)
thrbTxSeqs

writeXStringSet(thrbTxSeqs, "ENST00000396671.fasta")
