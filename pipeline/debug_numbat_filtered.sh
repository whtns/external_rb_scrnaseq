#!/bin/bash

/usr/bin/R CMD BATCH --no-restore --no-save \
"--args seu_path='/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat/SRR13884242_filtered_seu.rds' \
ref_path='/dataVolume/storage/Homo_sapiens/numbat/sridhar_ref.rds' \
t='1e-5' gamma='50' tau='0.3' init_k='6' read_prop='1' max_iter='2' min_LLR='2' cell_ceiling='1e4' \
max_entropy='0.8' allele_df='/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/numbat/SRR13884242_allele_counts.tsv.gz' \
matrix_file='/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/cellranger/SRR13884242/outs/filtered_feature_bc_matrix/matrix.mtx.gz' \
out_dir='/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/numbat_sridhar//SRR13884242' \
ncores='4' rprof_out='/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/numbat_sridhar/SRR13884242/log.prof'" \
scripts/run_numbat_filtered.R /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/logs/numbat_sridhar_SRR13884242.log
