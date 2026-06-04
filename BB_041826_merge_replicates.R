## I want to merge technical replicates for Kevin's dataset
## 4240, 4241 # no need to combine is uninformative
## 4242, 4243 # 2p balanced gain in one but not the other so nice to combine plus has a 1q+
## 4244, 4245 # no need incorrect 1q+ identified


## Merge 4240, 4241

# get cell ranger out put first
echo "sample_id,molecule_h5
SRR13884242,/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/cellranger/SRR13884242/outs/molecule_info.h5
SRR13884243,/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/cellranger/SRR13884243/outs/molecule_info.h5" > 4242_4243_aggr.csv

cellranger aggr \
--id=SRR13884242_SRR13884243_merged \
--csv=4242_4243_aggr.csv \
--normalize=mapped \
--localcores=8 \
--localmem=64



# For rep1: append "-1" to every CB tag
samtools view -h /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/cellranger/SRR13884242/outs/possorted_genome_bam.bam | \
awk 'BEGIN{OFS="\t"} /^@/{print; next} {
    for(i=12;i<=NF;i++) {
      if($i ~ /^CB:Z:/) $i = $i "-1"
    }
    print
  }' | \
samtools view -bS -o rep1_tagged.bam

# For rep2: append "-2" to every CB tag
samtools view -h /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/cellranger/SRR13884243/outs/possorted_genome_bam.bam | \
awk 'BEGIN{OFS="\t"} /^@/{print; next} {
    for(i=12;i<=NF;i++) {
      if($i ~ /^CB:Z:/) $i = $i "-2"
    }
    print
  }' | \
samtools view -bS -o rep2_tagged.bam

# Next merge BAM files
samtools merge -f SRR13884242_SRR13884243_merged.bam rep1_tagged.bam rep2_tagged.bam
# before running sort confirm if the bam file is sorted or not using
samtools view -H SRR13884242_SRR13884243_merged.bam | grep "@HD"
# if sorting is needed
samtools sort -o SRR13884240_SRR13884241_merged_sorted.bam SRR13884240_SRR13884241_merged.bam
samtools index SRR13884242_SRR13884243_merged.bam