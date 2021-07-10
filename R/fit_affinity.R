#' Fitting sigmoidal dose-response curves
#'
#' Fits a Kd curve for the equilibrium \eqn{[RL] <=> [R] + [L]},
#' where \eqn{[RL]} is measured by a proxy variable (response) in dependence of
#' the concentration \eqn{[L]} (dose).
#'
#' @param x A data frame to be evaluated.
#' @param formula A formula using variables of \code{x}; given as \code{RL ~ L0},
#' where \code{RL} is the variable that describes the (normalized) concentration
#' of the complex \eqn{[RL]_0}, and \code{L0} the concentration of titrated species
#' \eqn{[L]_0}.
#' @param R0 The total concentration of the species that is held constant during
#' the titration. If unknown (\code{NaN}), the constant ligand approximation is
#' used for fitting by default. See Details.
#' @param include_hill If \code{TRUE}, a Hill coefficient is applied onto \code{R};
#' not always appropiate.
#' @param limits_lower,limits_upper,limits_hill,limits_K_d Upper and lower bounds
#' for the paramters to be fitted.
#' @param ... Other arguments passed to the fitting algorithm \code{FUN}.
#' @param FUN The function used as fitting algorithm; must return a non-linear model.
#' Can be an unquoted function, but must be a character if prefixed as "package::function".
#' @param start A named list or named numeric vector of starting estimates. Guessed
#' from the data if \code{"auto"}; not passed if \code{NULL}.
#'
#' @details
#' If \code{R0} is a numeric value, the total concentration of the species that
#' is not titrated in the assay, the following quadratic equation is fitted:
#'
#' \deqn{[RL] = 1/(2 [R]_0) * ([R]_0 + [L]_0 + K_d - (([R]_0 + [L]_0 + K_d)^2 - 4 [R]_0 [L]_0)^(1/2))}
#'
#' else
#'
#' \deqn{[RL] = ([R]_0 [L]_0) / (K_d + [L]_0)}
#'
#' The latter makes use of the approximation \eqn{[L]_0 \approx [L]}, which is
#' valid only if \eqn{[R]_0 << K_d} over the inflection point of the curve. Else,
#' the data is biased by ligand depletion.
#'
#' In both cases, the fit can be normalized on the \eqn{[RL]} axis using
#' the fitted parameters \code{lower} and \code{upper} as \eqn{lower + (upper - lower) * [RL]}.
#'
#' The choice of the designations \eqn{R} (as mnemonic for "receptor") and
#' \eqn{L} (for "ligand") are arbitrary; the fitting works also in reversed
#' scenarios.
#'
#' \subsection{Choice of limits}{
#'
#' By default, the search space for the parameters, except for \code{K_d}, is
#' unconstraint (\code{c(-Inf, +Inf)}). This is because the fitting algorithm
#' may need to pass via "unreasonable" values during its search for the local/global
#' minimum.
#'
#' }
#'
#' \subsection{Choice of the fitting algorithm}{
#'
#' By default, \code{\link[minpack.lm:nlsLM]{minpack.lm::nlsLM}} is used as an
#' alternative to base R's \code{\link[stats:nls]{nls}}. You may consider using
#' \code{nls.multstart::nls_multstart} with additional arguments via \code{...}
#' to probe the dependence on different starting parameters. (See Examples.)
#'
#' }
#'
#' \subsection{Choice of starting values}{
#'
#' Starting values are guessed from the data if \code{"auto"}; by default:
#' \code{lower} is zero,
#' \code{upper} is the maximum of \eqn{[RL]} (LHS of \code{formula}),
#' \code{hill}, if present, is 1,
#' \code{K_d} is value of \eqn{[L]_0} (RHS of \code{formula})
#' which is closest to the half-maximal value of \eqn{[RL]}.
#'
#' You must set \code{start = NULL} with certain \code{FUN}, e.g. when using
#' \code{FUN = "nls.multstart::nls_multstart"}.
#'
#' }
#'
#' @examples
#' library(dplyr)
#' assay_file <- system.file("extdata", "gel_01.txt", package = "summerrband")
#' assay_data <- iqtl_meta(assay_file, list(conc = c(2^seq(10, 0), 0)))
#'
#' assay_data %>%
#'   filter(band_id == "band_1") %>%
#'   fit_Kd(vol_frac ~ conc, include_hill = FALSE, R0 = 2) %>%
#'   broom::tidy()
#'
#' assay_data %>%
#'   group_by(band_id) %>%
#'   summerr::model_cleanly_groupwise(fit_Kd, formula = vol_frac ~ conc)
#'
#' library(ggplot2)
#'
#' assay_data %>%
#'   group_by(band_id) %>%
#'   summerr::model_cleanly_groupwise(fit_Kd, formula = vol_frac ~ conc) %>%
#'   summerr::model_display(color = band_id) +
#'   scale_x_log10()
#'
#' # Other fitting algorithms
#'
#' assay_data %>%
#'   filter(band_id == "band_1") %>%
#'   fit_Kd(vol_frac ~ conc, include_hill = FALSE, R0 = 2, FUN = "nls", algorithm = "port") %>%
#'   broom::tidy()
#'
#' assay_data %>%
#'   filter(band_id == "band_1") %>%
#'   fit_Kd(vol_frac ~ conc, include_hill = FALSE, R0 = 2,
#'     FUN = "nls.multstart::nls_multstart", start = NULL,  ## important!
#'     iter = 500,
#'     start_lower = c(lower = 0,   upper = 0.5, K_d = 1e-3),
#'     start_upper = c(lower = 0.5, upper = 1.0, K_d = 1e3)) %>%
#'   broom::tidy()
#'
#' @export
fit_Kd <- function(x, formula, R0 = NaN, include_hill = FALSE,
                   limits_lower = c(-Inf, +Inf), limits_upper = c(-Inf, +Inf),
                   limits_hill = c(-Inf, +Inf), limits_K_d = c(0, 1e3),
                   start = "auto", ..., FUN = "minpack.lm::nlsLM") {

  RL <- formula.tools::lhs(formula)
  L0 <- formula.tools::rhs(formula)

  if (is.finite(R0)) {

    FML <- substitute(RL ~ I(lower + (upper - lower) * ((R0 + L0^hill + K_d) - sqrt(
      (R0 + L0^hill + K_d)^2 - 4 * R0 * L0^hill)) / (2 * R0)),
      list(R0 = R0, RL = RL, L0 = L0))

  } else {

    FML <- substitute(RL ~ I(lower + (upper - lower) * (L0^hill / (K_d + L0^hill))),
                      list(RL = RL, L0 = L0))

  }

  if (!is.null(start) && start == "auto") {

    starts <- list(hill = 1, lower = 0, upper = max(x[[RL]],
                                                    na.rm = TRUE),
                   K_d = x[[L0]][which.min(abs(
                     x[[RL]] - 0.5 * max(x[[RL]], na.rm = TRUE)))])

  }

  params <- list(hill = limits_hill, lower = limits_lower, upper = limits_upper,
                 K_d = limits_K_d)

  if (!include_hill) {

    starts$hill <- params$hill <- NULL
    FML <- do.call("substitute", list(FML, list(hill = 1)))

  }

  ll <- sapply(params, min, na.rm = TRUE, USE.NAMES = TRUE)
  ul <- sapply(params, max, na.rm = TRUE, USE.NAMES = TRUE)

  if (is.character(FUN)) {

    pts <- rev(strsplit(FUN, ":::?")[[1]])

    FUN <- pts[1]
    PKG <- pts[2]

    if (is.na(PKG)) rm(PKG)

  }

  if (is.null(start)) {

    # this allows compatibility with "nls.multstart::nls_multstart"

    eval(rlang::call2(FUN, x, formula = stats::as.formula(FML),
                      lower = ll, upper = ul,
                      ..., .ns = get0("PKG")))

  } else {

    eval(rlang::call2(FUN, x, formula = stats::as.formula(FML),
                      lower = ll, upper = ul,
                      start = starts, ..., .ns = get0("PKG")))

  }

}
