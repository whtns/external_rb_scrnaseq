
# All R scripts in ./R/ are sourced below, including packages.R and functions.R if present.
## Load your packages, e.g. library(targets).
suppressPackageStartupMessages(source("./packages.R"))

# --- Pipeline config flags (read by pipeline_targets_*.R files on source) ---
# Set to TRUE to stop tracking numbat RDS file content; prevents downstream
# rebuilds when numbat reruns update file timestamps/hashes without changing
# the sample set. Flip back to FALSE when you want the pipeline to detect
# genuinely new numbat outputs.
freeze_rds_files <- FALSE

## Load pipeline definition and constant files (functions are loaded via library(numbatHelpers))
lapply(list.files("./R", pattern = "^(pipeline_|constants)", full.names = TRUE), source)

# Seurat uses future internally for parallel ops; large Seurat objects exceed the
# 500 MB default limit. Remove the cap so workers can process them.
options(future.globals.maxSize = Inf)

# Absolute user library path — used both in script_lines (for R_LIBS_USER) and
# in tar_option_set(lib.loc=) so workers find packages regardless of HOME resolution.
.user_lib <- "/home1/stachele/R/x86_64-pc-linux-gnu-library/4.4"

# Lines added to every SLURM worker script: load modules and set library paths.
# Workers start as fresh R sessions on compute nodes, so they need this environment.
.slurm_script_lines <- c(
  "source /etc/profile.d/modules.sh",
  "module load r/4.4.1 curl bzip2 libxml2 cairo fontconfig freetype poppler/23.04.0",
  paste0("export LD_LIBRARY_PATH=/home1/stachele/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"),
  paste0("export R_LIBS_USER=", .user_lib)
)

# Light controller: file-tracking, plotting, metadata (~8 GB total, 2 CPUs).
.worker_light <- crew_controller_slurm(
  name                        = "light",
  workers                     = 8,
  seconds_idle                = 120,
  seconds_timeout             = 600,
  slurm_partition             = "epyc-64",
  slurm_cpus_per_task         = 2,
  slurm_memory_gigabytes_per_cpu = 4,
  slurm_time_minutes          = 240,
  slurm_log_output            = "logs/crew_light_%j.out",
  slurm_log_error             = "logs/crew_light_%j.err",
  script_lines                = .slurm_script_lines
)

# Heavy controller: Seurat processing, numbat, integration, diffex (~64 GB total, 4 CPUs).
.worker_heavy <- crew_controller_slurm(
  name                        = "heavy",
  workers                     = 4,
  seconds_idle                = 300,
  seconds_timeout             = 600,
  slurm_partition             = "epyc-64",
  slurm_cpus_per_task         = 4,
  slurm_memory_gigabytes_per_cpu = 16,
  slurm_time_minutes          = 720,
  slurm_log_output            = "logs/crew_heavy_%j.out",
  slurm_log_error             = "logs/crew_heavy_%j.err",
  script_lines                = .slurm_script_lines
)

tar_option_set(
  memory             = "transient",
  garbage_collection = TRUE,
  error              = "continue",
  workspace_on_error = TRUE,
  trust_timestamps   = TRUE,
  storage            = "main",
  library            = .user_lib,
  controller         = crew_controller_group(.worker_heavy, .worker_light),
  # Default to heavy; workers are persistent so lightweight targets on a heavy worker
  # cost nothing extra. Tag individual targets with controller = "light" to opt down.
  resources          = tar_resources(crew = tar_resources_crew(controller = "heavy"))
)

## _targets.R must return a list of tar_target objects.
## Each pipeline_targets_* variable is a list of targets; c() flattens them.
c(
  pipeline_targets_inputs,
  pipeline_targets_seurat,
  pipeline_targets_diffex,
  pipeline_targets_integration,
  pipeline_targets_figures,
  # Shallow diagnostics on the numbat RDS inputs. Depends only on
  # numbat_rds_srx, so it can be built alone without waking the rest.
  pipeline_targets_qc
)

# Sample-level analysis notes for the SRR era moved to
# docs/srr_sample_analysis_notes_2023.md (comments only; nothing read them).
