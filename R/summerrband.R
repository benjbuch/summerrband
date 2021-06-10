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
#' \subsection{Importing IQTL Data}{
#'
#' \link{iqtl_read}, \link{iqtl_meta}, \link{iqtl_import_all}
#'
#' \link{iqtl_view} to quickly view an IQTL file
#'
#' }
#'
#' \subsection{Workflows}{
#'
#' See \link{use_template} for workflows involving, e.g., fitting of Kd values (\link{fit_Kd}).
#'
#' }
#'
#' @details
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
#'
#' \item{A01}{EMSA with selectivity measurements}
#'
#' \item{A02}{EMSA with affinity measurements and Kd fitting}
#'
#' }
#'
#' @export
use_template <- function(version = "A01") {

  summerr::get_template(package = "summerrband", filename = "template",
                        version = version)

}
