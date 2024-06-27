library(biomaRt)

mart <- useMart("ensembl",
                dataset="hsapiens_gene_ensembl")

seq = getSequence(id = "ENST00000396671",
                  type = "ensembl_transcript_id",
                  seqType = "cdna",
                  mart = mart)
show(seq)
