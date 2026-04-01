Summary of changes:

pipeline/config.yaml:149 — added cell_ceiling_by_sample: SRR27187901: 1e4
pipeline/Snakefile:139 — added cell_ceiling_for_sample() helper (mirrors subset_bad_cell_types_for_sample pattern)
pipeline/Snakefile:1327 — numbat_sridhar_filtered now calls the helper instead of the global config
To add overrides for other samples in the future, just add entries under cell_ceiling_by_sample in config.yaml.