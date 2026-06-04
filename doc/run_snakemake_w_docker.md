# Run a merged sample
docker run --rm \
  --memory=55g --memory-swap=55g \
  --ulimit nofile=65536:65536 \
  -v /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj:/proj \
  -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
  -w /proj/pipeline \
  numbat-pipeline \
  snakemake --snakefile Snakefile --cores 6 --rerun-incomplete \
  --config run_samples=SRR13884242_SRR13884243_merged


## Unlock after interrupted run (run this first if Snakemake reports a lock)
docker run --rm \
  --memory=55g --memory-swap=55g \
  --ulimit nofile=65536:65536 \
  -v /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj:/proj \
  -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
  -w /proj/pipeline \
  numbat-pipeline \
  snakemake --snakefile Snakefile --unlock

## Run a single sample (pass run_samples via --config)
docker run --rm \
  --memory=55g --memory-swap=55g \
  --ulimit nofile=65536:65536 \
  -v /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj:/proj \
  -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
  -w /proj/pipeline \
  numbat-pipeline \
  snakemake --snakefile Snakefile --cores 6 --rerun-incomplete --config run_samples=SRR17960482 \
  --cores 6 --dry-run

## Dry-run (check what would execute)
docker run --rm \
  --memory=55g --memory-swap=55g \
  --ulimit nofile=65536:65536 \
  -v /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj:/proj \
  -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
  -w /proj/pipeline \
  numbat-pipeline \
  snakemake --snakefile Snakefile --cores 6 --dry-run

## Full run
docker run --rm \
  --memory=55g --memory-swap=55g \
  --ulimit nofile=65536:65536 \
  -v /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj:/proj \
  -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
  -w /proj/pipeline \
  numbat-pipeline \
  snakemake --snakefile Snakefile --cores 6

## Force re-run pileup (use after fixing SNP reference chr naming)
docker run --rm \
  --memory=55g --memory-swap=55g \
  --ulimit nofile=65536:65536 \
  -v /dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj:/proj \
  -v /dataVolume/storage/Homo_sapiens:/dataVolume/storage/Homo_sapiens:ro \
  -w /proj/pipeline \
  numbat-pipeline \
  snakemake --snakefile Snakefile --cores 6 --forcerun pileup
