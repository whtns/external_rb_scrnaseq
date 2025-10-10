check_tcga_peaks <- function(cohort, arm_ranges = list(
  "1q" = list(chrom = "1", start = 121700000, end = 248956422, direction = "Amp"),
  "2p" = list(chrom = "2", start = 1, end = 93400000, direction = "Amp"),
  "6p" = list(chrom = "6", start = 1, end = 61000000, direction = "Amp"),
  "16q" = list(chrom = "16", start = 33100000, end = 90354753, direction = "Del")
)) {

    # browser()

    tcga_cohort <- 
    cohort |>
    TCGAgistic::tcga_gistic_load(source = "Firehose", cnLevel = "all")


    df <- tcga_cohort@cytoband.summary |>
    separate(Wide_Peak_Limits, into = c("chrom", "range"), sep = ":") %>%
    separate(range, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.integer(start), end = as.integer(end)) |>
    mutate(chrom = str_remove(chrom, "chr"))
    found <- sapply(names(arm_ranges), function(arm) {
      arm_info <- arm_ranges[[arm]]
    #   browser()
      any(
        df$chrom == arm_info$chrom &
        df$qvalues < 0.05 &
        df$start <= arm_info$end &
        df$end >= arm_info$start &
        df$Variant_Classification == arm_info$direction
      )
    })
    data.frame(abbreviation = cohort, t(found))
}

check_tcga_peaks(tcga_cohorts[[1]])

test1 <- 
tcga_cohorts |> 
  map(check_tcga_peaks) |>
  dplyr::bind_rows(.id = "cohort") |>
  identity()
