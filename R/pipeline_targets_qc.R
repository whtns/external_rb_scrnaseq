# Diagnostic targets for the numbat RDS inputs.
# Defines: pipeline_targets_qc (spliced into the target list in _targets.R)
#
# These targets exist to answer one question cheaply: is each *_numbat.rds
# object what we think it is? scripts/process_numbat_rds.R now builds objects at
# a per-sample consensus round (results/numbat_selected_round.csv) rather than
# numbat's package default i = 2, and nothing about the file's name or location
# records which round it got. So this reads each RDS, reports the round it
# actually holds, and renders numbat's own two plots for eyeballing.
#
# Deliberately shallow. It depends only on numbat_rds_srx -- no Seurat objects,
# no cluster dictionary, no clone simplifications -- so it can be built without
# touching the expensive half of the pipeline, and so a failure here implicates
# the numbat inputs rather than anything downstream.
#
# Build it on its own:
#   tar_make(names = c(numbat_rds_qc_table, numbat_rds_qc_report))
# (see qc_numbat_rds.sbatch). Everything is read-only with respect to
# output/numbat_sridhar/; the only writes are under results/numbat_round_qc/.

pipeline_targets_qc <- list(

  # SRX only. The SRR samples are out of scope for this project (AGENTS.md:
  # never rebuild them), nothing else in the pipeline maps over numbat_rds_all,
  # and including them just buries the samples we care about in a 71-row table.
  # Their objects are untouched on disk; they are simply not inspected here.
  tarchetypes::tar_files(
    numbat_rds_srx,
    {
      f <- retrieve_numbat_rds_files("output/numbat_sridhar/")
      f[grepl("^SRX", basename(f))]
    },
    format = "file"
  ),

  # One branch per RDS. Each branch opens its file exactly once, writes that
  # sample's QC PDF as a side effect, and returns a one-row summary carrying the
  # PDF path -- cheaper than a separate plotting target, which would re-read
  # objects that run to 460 MB on disk.
  #
  # Heavy controller on purpose: bulk_clones alone is millions of rows on the
  # larger samples, so the light controller's 8 GB is not enough headroom.
  tar_target(
    numbat_rds_qc,
    numbatHelpers::qc_numbat_rds(
      numbat_rds_srx,
      out_dir  = "results/numbat_round_qc",
      manifest = "results/numbat_selected_round.csv"
    ),
    pattern   = map(numbat_rds_srx),
    iteration = "vector",
    resources = .heavy_resources
  ),

  # Collated table. Written as a file target so it can be opened outside the
  # pipeline and diffed between runs.
  tar_target(
    numbat_rds_qc_table,
    {
      dir.create("results", showWarnings = FALSE)
      qc <- numbat_rds_qc[order(numbat_rds_qc$sample_id), ]
      out <- "results/numbat_rds_qc.csv"
      readr::write_csv(qc, out)
      out
    },
    format    = "file",
    resources = .light_resources
  ),

  # Every per-sample PDF in one file, ordered by sample, for flipping through.
  # Branches that failed to render contribute no page rather than a blank one.
  tar_target(
    numbat_rds_qc_report,
    {
      qc  <- numbat_rds_qc[order(numbat_rds_qc$sample_id), ]
      pdfs <- qc$pdf[!is.na(qc$pdf) & file.exists(qc$pdf)]
      out <- "results/numbat_rds_qc_report.pdf"
      if (length(pdfs) == 0) {
        # Still produce the file so the target has a stable output; an empty
        # report is itself the diagnostic.
        grDevices::pdf(out, width = 8, height = 4)
        plot.new(); title("no numbat QC pages rendered")
        grDevices::dev.off()
      } else {
        qpdf::pdf_combine(pdfs, output = out)
      }
      out
    },
    format    = "file",
    resources = .light_resources
  )
)
