NULL

#' Separate TCGA barcode
#'
#' A wrapper around \link[tidyr]{separate} to split TCGA barcodes into individual variables as described \href{https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#reading-barcodes}{here}.
#'
#' @importFrom dplyr select
#' @importFrom tidyr separate
#' @importFrom stringr str_count str_detect
#' @importFrom rlang enquo quo_name
#'
#' @param data A data frame.
#' @param col Column name containing the barcode
#' @param remove If TRUE, remove the column containing the barcode.
#'
#' @export
separate_barcode <- function(data, col = barcode, remove = TRUE) {
  col <- quo_name(enquo(col))
  barcode <- data[[col]]

  n_elts <- unique(str_count(barcode, "-")) + 1
  if (length(n_elts) != 1) stop("all barcodes should be of same length")

  barcode_regex <- c("^TCGA", "\\w{2}", "\\w{4}", "\\d{2}[A-Z]", "\\d{2}[A-Z]?", "\\w{4}", "\\d{2}")
  barcode_regex <- paste0(paste(barcode_regex[1:n_elts], collapse = "-"), "$")
  valid_barcode <- str_detect(barcode, barcode_regex)
  if (!all(valid_barcode)) stop("invalid barcode found: ", glue_collapse(barcode[!valid_barcode], sep = ", "))

  into <- c(".project", "tss", "participant", "sample_vial", "portion_analyte", "plate", "center")

  if (any(into %in% colnames(data))) stop(colnames(data)[colnames(data) %in% into], " column already exists")

  data <- separate(data, col, into[1:n_elts], sep = "-", remove = remove, fill = "right")
  if (n_elts > 3) data <- separate(data, "sample_vial", c("sample", "vial"), 2)
  if (n_elts > 4) data <- separate(data, "portion_analyte", c("portion", "analyte"), 2)
  data <- mutate_at(data, "analyte", str_replace, "^$", NA_character_)
  select(data, -".project")
}

#' Get the expression values out of the TCGA SummarizedExperiment object
#'
#' A wrapper around \link[SummarizedExperiment]{assay} to extract the TCGA assay, filter out the genes of interest and
#' always return gene_symbols instead of Ensembl gene IDs.
#'
#' @param data A SummarizedExperiment object.
#' @param goi A character vector containing the gene symbols of interest.
#' @param ... Extra arguments
#'
#' @export
get_tcga_goi <- function (data, goi, ...) {
  UseMethod("get_tcga_goi", data)
}


#' @rdname get_tcga_goi
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr gather
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @importFrom geneid ensgid_to_gene_symbol
#'
#' @method get_tcga_goi RangedSummarizedExperiment
#' @export
get_tcga_goi.RangedSummarizedExperiment <- function(data, goi, ...) {
  x <- assay(data)
  x <- as_tibble(x, rownames = "gene_symbol")
  if (all(str_detect(x[["gene_symbol"]], "^ENSG[0-9]{11}$"))) {
    x["gene_symbol"] <- ensgid_to_gene_symbol(x[["gene_symbol"]])
  }
  x <- filter(x, gene_symbol %in% goi)
  x <- gather(x, "barcode", "expression", -gene_symbol)
  x
}

#' @rdname get_tcga_goi
#'
#' @param element Index of name of the list element to use.
#'
#' @method get_tcga_goi list
#' @export
get_tcga_goi.list <- function(data, goi, element = 1, ...) {
  if (missing(element) && length(data) > 1) warning("Using first element of the list")
  get_tcga_goi(data[[element]], goi)
}

#' @rdname get_tcga_goi
#'
#' @method get_tcga_goi character
#' @export
get_tcga_goi.character <- function(data, goi, ...) {
  get_tcga_goi(read_rds(data), goi, ...)
}
