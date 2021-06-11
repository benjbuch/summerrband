#' Import an IQTL file
#'
#' Reads a file exported with ImageQuant TL and performs basic actions to tidy
#' up the format.
#'
#' @param file Path to the file.
#' @param row_name Name to give to rows.
#' @param col_group Name to give to column groups, i.e., repeating units of \code{col_repeat}.
#' @param col_repeat Column variables that repeat along \code{col_group}.
#' @param .id General identifier added to \code{row_name} and \code{col_group}.
#' @param .sep General separator.
#' @param reverse_rows If \code{TRUE}, the "\code{row_name}" identifiers are counted
#' from bottom to top; else from top to bottom.
#' @param gel_name A \link[base:regex]{regular expression} to extract an identifier
#' from the file name.
#'
#' @return
#' A data frame of class \code{\link[data.table:data.table]{data.table}}.
#'
#' @examples
#' assay_file <- system.file("extdata", "gel_01.txt", package = "summerrband")
#'
#' iqtl_read(assay_file)
#'
#' @importFrom data.table :=
#'
#' @export
iqtl_read <- function(file, row_name = "band", col_group = "lane",
                      col_repeat = c(quantify = "vol", "vol+bkg", "CAL", "MW", "Rf"),
                      reverse_rows = TRUE, .id = "id", .sep = "_",
                      gel_name = paste0("gel", .sep, "[[:alnum:]_-]+")) {

  # visible binding for data.table variables
  .SD <- NULL

  # number of columns per lane in IQTL file

  rows = c(paste(col_group, .id, sep = .sep))
  cols = c(paste(row_name,  .id, sep = .sep), col_repeat)

  # check if decimal comma or decimal point was used in column to quantify ...

  decimal_point <- sapply(data.table::fread(file, skip = 1, header = TRUE,
                                            na.strings = "-", nrows = 1, dec = "."),
                          is.numeric)[which(names(cols) == "quantify")]

  # import IQTL file as is; skip the first line since it contains but the lane
  # number information that can be inferred

  data_raw <- data.table::fread(file, skip = 1, header = TRUE, na.strings = "-",
                                dec = c(",", ".")[decimal_point + 1])

  data_tmp <- data.table::data.table(matrix(t(as.matrix(data_raw)), ncol = length(cols),
                                            byrow = TRUE, dimnames = list(NULL, cols)))

  rowvalues <- data_tmp[[paste(row_name, .id, sep = .sep)]]

  if (reverse_rows) {

    # assign band_id from bottom to top, which is more convenient to work with

    data.table::set(data_tmp, i = NULL, j =  paste(row_name, .id, sep = .sep),
                    value = paste(row_name, max(rowvalues) - rowvalues, sep = .sep))

  } else {

    # assign band_id_iqtl from top to bottom, which is the default in IQTL

    data.table::set(data_tmp, i = NULL, j =  paste(row_name, .id, sep = .sep),
                    value = paste(row_name, rowvalues, sep = .sep))

  }

  colname <- paste(col_group, seq(dim(data_raw)[[2]] / length(cols)), sep = .sep)

  data.table::set(data_tmp, i = NULL, j = rows, value = factor(
    rep(colname, dim(data_raw)[[1]]), levels = colname))

  # keep file name and extract gel identifier if possible

  data.table::set(data_tmp, i = NULL, j = "file", value = file)

  data.table::set(data_tmp, i = NULL, j = paste0("gel", .sep, .id),
                  value = stringr::str_extract(file, gel_name))

  # calculate relative lane composition

  data_tmp[, paste0(cols[["quantify"]], .sep, "frac") := prop.table(.SD),
           by = c(paste0(col_group, .sep, .id)),
           .SDcols = cols[["quantify"]]]

  data.table::setcolorder(
    data_tmp,
    neworder = c(rows[[1]], cols[[1]], cols[["quantify"]],
                 paste0(cols[["quantify"]], .sep, "frac"),
                 setdiff(cols, c(cols[[1]], cols[["quantify"]],
                                 paste0(cols[["quantify"]], .sep, "frac"))),
                 "file"))

  data.table::setkeyv(data_tmp, c(rows[[1]], cols[[1]]))

  data_tmp[]  # guarantee to print

}

#' Import and view an IQTL file
#'
#' Reads a file exported with ImageQuant TL and performs basic actions to tidy
#' up the format using \code{\link{iqtl_read}} and draws an overview plot.
#'
#' @inheritParams iqtl_read
#' @param ... Other arguments passed to \code{iqtl_read}.
#' @param selected_bands Character vector to imit plotting to specific bands,
#' e.g. \code{"band_0"}.
#' @param on_y_axis Which variable to draw on the y axis.
#'
#' @details
#' This function calls \code{iqtl_read} with its default arguments for \code{row_name}
#' ("band"), \code{col_group} ("lane"), \code{.id} ("id") and \code{.sep} ("_").
#' It returns the same object as a default call to \code{iqtl_read}.
#' As a side-effect, a plot is printed to the current device.
#'
#' @return
#' A data frame of class \code{\link[data.table:data.table]{data.table}}.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
iqtl_view <- function(file, selected_bands = character(), ..., on_y_axis = "vol_frac") {

  gel <- iqtl_read(file = file, row_name = "band", col_group = "lane", .id = "id",
                   .sep = "_", ...)

  on_y_axis <- rlang::ensym(on_y_axis)

  plt <- gel %>%
    dplyr::filter(.data$band_id %in% if (length(!!selected_bands) == 0) .$band_id else !!selected_bands) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$lane_id, y = !!on_y_axis,
                                 shape = .data$band_id)) +
    ggplot2::geom_point() +
    ggplot2::geom_path(ggplot2::aes(group = .data$band_id)) +
    ggplot2::labs(title = unique(gel$gel_id), caption = file) +
    ggplot2::theme_linedraw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  print(plt)

  gel[]  # guarantee to print

}

#' Associate metadata with an IQTL file
#'
#' Add metadata along an imported IQTL data file.
#'
#' @inheritParams iqtl_read
#' @param ... Other arguments passed to \code{iqtl_read}.
#' @param meta_data List to construct a metadata data frame from.
#' @param exclude A character vector with \code{lane_id}s to remove from the result.
#'
#' @details
#' This function calls \code{iqtl_read} with its default arguments for
#' \code{col_group} ("lane"), \code{.id} ("id") and \code{.sep} ("_").
#'
#' If \code{meta_data} contains an entry named \code{"lane_id"}, its entries are
#' used to join the imported data with the metadata annotation. If no such entry
#' exists, the data is assumed to be parallel to the result of \code{iqtl_read}.
#' (More specifically, to the \code{levels} of \code{.$lane_id}.)
#'
#' Internally, \code{meta_data} is converted to a \code{data.frame}.
#'
#' @return
#' A data frame of class \code{\link[data.table:data.table]{data.table}}.
#'
#' @examples
#' assay_file <- system.file("extdata", "gel_01.txt", package = "summerrband")
#' assay_data <- iqtl_meta(assay_file, list(conc = c(2^seq(10, 0), 0)))
#'
#' @export
iqtl_meta <- function(file, meta_data = list(conc = NA), exclude = NULL, ...) {

  gel <- iqtl_read(file = file, col_group = "lane", .id = "id", .sep = "_", ...)

  if ("lane_id" %in% names(meta_data)){

    y <- dplyr::left_join(gel, data.frame(meta_data), by = "lane_id")

  } else {

    y <- dplyr::left_join(gel, data.frame(meta_data,
                                          lane_id = levels(gel$lane_id)),
                          by = "lane_id")

  }

  # exclude lanes when "exclude" given

  y[!(y$lane_id %in% exclude), ][]

}

#' Import IQTL files from a list
#'
#' Given a \code{list} or \code{data.frame} with an entry/column \code{file},
#' this function feeds its arguments to \code{\link{iqtl_meta}}.
#'
#' @param l A list with descriptors. See examples.
#' @param path A path in which all \code{file}s are residing. If \code{NULL},
#' \code{file}s should be fully expanded paths or refer to the current working
#' directory.
#' @inheritDotParams iqtl_meta
#'
#' @details
#' Entries named \code{.$file} or \code{.$exclude} are passed as separate arguments
#' to \code{\link{iqtl_meta}}. All other entries of \code{l} are passed to
#' \code{iqtl_meta} as if they were given as \code{meta_data} list.
#'
#' Empty entries are dropped without a warning. Missing files are dropped with
#' a warning.
#'
#' @examples
#'
#' @export
iqtl_import_all <- function(l, path = NULL, ...) {

  # visible binding for purrr variables
  .x <- NULL

  l <- purrr::compact(l)  # drop empty entries

  if (!is.null(path)) l <- purrr::map(l, ~ purrr::modify_in(
    .x, "file", ~ file.path(path, .x)))

  purrr::map_dfr(l, ~ iqtl_meta(file = .x$file, meta_data = .x[setdiff(
    names(.x), c("file", "exclude"))]), exclude = .x$exclude, ... = ...)

}
