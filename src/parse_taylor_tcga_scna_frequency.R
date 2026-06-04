library(tidyverse)
library(fs)
library(readxl)
library(janitor)

#' Reorder an x or y axis within facets
#'
#' Reorder a column before plotting with faceting, such that the values are ordered
#' within each facet. This requires two functions: \code{reorder_within} applied to
#' the column, then either \code{scale_x_reordered} or \code{scale_y_reordered} added
#' to the plot.
#' This is implemented as a bit of a hack: it appends ___ and then the facet
#' at the end of each string.
#'
#' @param x Vector to reorder.
#' @param by Vector of the same length, to use for reordering.
#' @param within Vector of the same length that will later be used for faceting
#' @param fun Function to perform within each subset to determine the resulting
#' ordering. By default, mean.
#' @param sep Separator to distinguish the two. You may want to set this manually
#' if ___ can exist within one of your labels.
#' @param ... In \code{reorder_within} arguments passed on to \code{\link{reorder}}.
#' In the scale functions, extra arguments passed on to
#' \code{\link[ggplot2]{scale_x_discrete}} or \code{\link[ggplot2]{scale_y_discrete}}.
#'
#' @source "Ordering categories within ggplot2 Facets" by Tyler Rinker:
#' \url{https://trinkerrstuff.wordpress.com/2016/12/23/ordering-categories-within-ggplot2-facets/}
#'
#' @examples
#'
#' library(tidyr)
#' library(ggplot2)
#'
#' iris_gathered <- gather(iris, metric, value, -Species)
#'
#' # reordering doesn't work within each facet (see Sepal.Width):
#' ggplot(iris_gathered, aes(reorder(Species, value), value)) +
#'   geom_boxplot() +
#'   facet_wrap(~ metric)
#'
#' # reorder_within and scale_x_reordered work.
#' # (Note that you need to set scales = "free_x" in the facet)
#' ggplot(iris_gathered, aes(reorder_within(Species, value, metric), value)) +
#'   geom_boxplot() +
#'   scale_x_reordered() +
#'   facet_wrap(~ metric, scales = "free_x")
#'
#' @export
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
	new_x <- paste(x, within, sep = sep)
	stats::reorder(new_x, by, FUN = fun)
}


#' @rdname reorder_within
#' @export
scale_x_reordered <- function(..., sep = "___") {
	reg <- paste0(sep, ".+$")
	ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}


#' @rdname reorder_within
#' @export
scale_y_reordered <- function(..., sep = "___") {
	reg <- paste0(sep, ".+$")
	ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

tcga_abbreviations <- 
	tibble::tribble(
		~Study.Abbreviation,                                                        ~Study.Name,
		"LAML",                                           "Acute Myeloid Leukemia",
		"ACC",                                         "Adrenocortical carcinoma",
		"BLCA",                                     "Bladder Urothelial Carcinoma",
		"LGG",                                         "Brain Lower Grade Glioma",
		"BRCA",                                        "Breast invasive carcinoma",
		"CESC", "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
		"CHOL",                                               "Cholangiocarcinoma",
		"LCML",                                     "Chronic Myelogenous Leukemia",
		"COAD",                                             "Colon adenocarcinoma",
		"CNTL",                                                         "Controls",
		"ESCA",                                             "Esophageal carcinoma",
		"FPPP",                                              "FFPE Pilot Phase II",
		"GBM",                                          "Glioblastoma multiforme",
		"HNSC",                            "Head and Neck squamous cell carcinoma",
		"KICH",                                               "Kidney Chromophobe",
		"KIRC",                                "Kidney renal clear cell carcinoma",
		"KIRP",                            "Kidney renal papillary cell carcinoma",
		"LIHC",                                   "Liver hepatocellular carcinoma",
		"LUAD",                                              "Lung adenocarcinoma",
		"LUSC",                                     "Lung squamous cell carcinoma",
		"DLBC",                  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
		"MESO",                                                     "Mesothelioma",
		"MISC",                                                    "Miscellaneous",
		"OV",                                "Ovarian serous cystadenocarcinoma",
		"PAAD",                                        "Pancreatic adenocarcinoma",
		"PCPG",                               "Pheochromocytoma and Paraganglioma",
		"PRAD",                                          "Prostate adenocarcinoma",
		"READ",                                            "Rectum adenocarcinoma",
		"SARC",                                                          "Sarcoma",
		"SKCM",                                          "Skin Cutaneous Melanoma",
		"STAD",                                           "Stomach adenocarcinoma",
		"TGCT",                                      "Testicular Germ Cell Tumors",
		"THYM",                                                          "Thymoma",
		"THCA",                                                "Thyroid carcinoma",
		"UCS",                                           "Uterine Carcinosarcoma",
		"UCEC",                             "Uterine Corpus Endometrial Carcinoma",
		"UVM",                                                   "Uveal Melanoma"
	) |> 
	clean_names()


taylor_freq0 <- read_xlsx("data/taylor_et_al_2018_genomic_and_functional_approaches_to_understanding_cancer_aneuploidy/1-s2.0-S1535610818301119-mmc2.xlsx", skip = 1) |> 
	clean_names() |> 
	identity()

taylor_freq1 <- 
	taylor_freq0 |> 
	tidyr::pivot_longer(starts_with("x"), names_to = "arm", values_to = "change") |> 
	dplyr::mutate(arm = str_remove(arm, "x")) |> 
	dplyr::filter(arm %in% c("1q", "2p", "6p", "16q")) |> 
	dplyr::mutate(rb_scna = case_when(
		(arm == "1q" & change == 1) ~ 1,
		(arm == "2p" & change == 1) ~ 1,
		(arm == "6p" & change == 1) ~ 1,
		(arm == "16q" & change == -1) ~ -1,
		.default = 0
	)) |> 
	dplyr::group_by(type, arm) |> 
	dplyr::mutate(arm = case_when(
		arm == "1q" ~ "1q+",
		arm == "2p" ~ "2p+",
		arm == "6p" ~ "6p+",
		arm == "16q" ~ "16q-",
	)) |> 
	# dplyr::summarize(percent_affected = sum(rb_scna)) |>
	dplyr::summarize(percent_affected = abs(sum(rb_scna)/dplyr::n())) |>
	dplyr::mutate(arm = factor(arm, levels = c("1q+", "2p+", "6p+", "16q-"))) |> 
	dplyr::arrange(arm, desc(percent_affected)) |> 
	dplyr::left_join(tcga_abbreviations, by = c("type" = "study_abbreviation")) |> 
	identity()

ggplot(taylor_freq1, aes(x = reorder_within(type, percent_affected, arm), y = percent_affected, fill = type)) + 
	geom_col() + 
	scale_x_reordered() +
	facet_wrap(~arm, scales = "free_x") + 
	theme_minimal() + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				axis.title.x = element_blank())

#' ggplot(iris_gathered, aes(reorder_within(Species, value, metric), value)) +
#'   geom_boxplot() +
#'   scale_x_reordered() +
#'   facet_wrap(~ metric, scales = "free_x")
