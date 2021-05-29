utils::globalVariables(".")

#' summerrband-package
#'
#' @author Benjamin Buchmuller \email{benjamin.buchmuller@@tu-dortmund.de}
#'
#' @title summerrband
#'
#' @description
#'
#' Analyze gel-shift assays.
#'
#' @details
#'
#' This package provides functions to evaluate gel-shift assays, e.g. to determine
#' binding affinity.
#'
#' @docType package
#' @name summerrband-package
NULL

#' Use a template for a IQTL/EMSA script
#'
#' @param version A template version identifier.
#'
#' @details
#'
#' \describe{
#' \item{A01}{EMSA for selectivity}
#' \item{A02}{EMSA for affinity}
#' }
#'
#' @export
use_template <- function(version = "A01") {

  summerr::get_template(package = "summerrband", filename = "template",
                        version = version)

}
