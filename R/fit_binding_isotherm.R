#' Get macroscopic binding polynomial
#'
#' Returns an expression that represents the macroscopic (thermodynamic) binding
#' polynomial (the grand canonical partition function) for an idealized system
#' that binds up to i ligands.
#'
#' @param params A character vector of length i to model i binding states; the
#' associated equilibrium constants take these names.
#' @param degree The degree of the polynom; overwritten by \code{params}.
#' @param kname The common variable name for the equilibrium constants.
#' @param xname The name of the variable representing the bound/titrated ligand.
#'
#' @return An expression.
#'
#' @export
gpf_macro <- function(degree = 2, params = as.character(seq(degree)), xname = "x",
                      kname = "K_") {

  # 1 + K_1 * x + K_1 * K_2 * x^2 + K_1 * K_2 * K_3 * x^3

  parse(text = paste0("1 + ", paste0(sapply(seq_along(params), function(m) paste0(
    paste0(kname, params[1:m], collapse = " * "), " * ", xname, "^", m)), collapse = " + ")))

}

#' Get microscopic binding polynomial
#'
#' Returns an expression that represents the microscopic binding
#' polynomial for an idealized system with i (different) binding sites.
#'
#' @param params A character vector of length i to model i binding sites; the
#' associated equilibrium constants take these indices.
#' @inheritParams gpf_macro
#'
#' @details
#' As compared to \code{\link{gpf_macro}}, the returned expression has
#' more variables than constructed from \code{degree} or passed via \code{params}.
#'
#' @return An expression.
#'
#' @export
gpf_micro <- function(degree = 2, params = letters[seq(degree)], xname = "x", kname = "k_") {

  # 1 + (k_a + k_b + k_c) * x + (k_ab + k_bc + k_ac) * x^2 + ...

  parse(text = paste0("1 + ", paste0(sapply(seq_along(params), function(m) paste0(
    "(", paste0(kname, utils::combn(
      params, m, FUN = paste0, collapse = ""), collapse = " + "), ") * ", xname, "^", m)),
    collapse = " + ")))

}

#' Compute the bound fraction from the grand partition function
#'
#' @param x A vector of ligand concentrations.
#' @param binding_constants The values of the association equilibrium constant
#' given as a (named) vector. See Details.
#' @param type Either \code{"micro"} or \code{"macro"} (default).
#'
#' @details
#' If \code{type == "macro"}, the binding isotherm is computed for a macroscopic
#' process using \code{\link{gpf_macro}}. The entries are assumed to be in ascending
#' order, i.e., \code{binding_constants[2]} gives the apparent thermodynamic binding
#' constant for associating the second ligand.
#'
#' If \code{type == "micro"} the vector must contain the (named) microsocopic
#' binding constants compatible with \code{\link{gpf_micro}}. Note that in this
#' case, the entries are assumed to be in combinatorically increasing order,
#' i.e., "a", "b", "c", "ab", "ac", "bc", "abc", for a third degree polynom.
#' The entries can be named, but the names are not used for the assignment due
#' to the complexity of the problem.
#'
#' In either case, the degree is determined from the length of the binding constants
#' vector.
#'
#' @export
gpf_fraction_bound <- function(x, binding_constants, type = "macro") {

  if (type == "macro") {

    # macroscopic binding polynomial

    binding_constants <- sort(binding_constants[which(binding_constants > 0)],
                              decreasing = TRUE, na.last = NA)

    bd <- length(binding_constants)

    if (is.null(names(binding_constants))) {

      kname <- "K_"
      bn <- paste0(kname, as.character(seq_along(binding_constants)))

    } else {

      names(binding_constants)[which(names(binding_constants) == "x")] <- "..x"

      kname <- " "
      bn <-  names(binding_constants)

    }

  } else {

    # microscopic binding polynomial; assume that a user will not provide more
    # than 1000 parameters here ...

    bd <- which(sapply(1:10, function(i) sum(choose(i, 1:i))) == length(binding_constants))

    if (length(bd) == 0) stop("Some binding constants are missing (or superfluous).")

    kname <- "k_"
    bn <- unlist(sapply(seq(bd), function(m) paste0(kname, utils::combn(
      letters[1:bd], m, FUN = paste0, collapse = ""))))

  }

  names(binding_constants) <- bn

  gpf <- do.call(paste0("gpf_", type), list(params = sub(
    kname, "", names(binding_constants)[1:bd], fixed = TRUE),
    xname = "x", kname = kname))
  gpf_dx <- stats::D(gpf, "x")

  summerr::log_debugging(gpf)

  if (getOption("summerr.debug", FALSE)) summerr::log_object(binding_constants)

  gpf_env <- c(list(x = x), as.list(binding_constants))

  res <- eval(gpf_dx, envir = gpf_env) * x / eval(gpf, envir = gpf_env)

  attr(res, "degree") <- bd
  attr(res, "params") <- binding_constants

  res

}

#' Quick overview plots
#'
#' @noRd
gpf_fraction_plot <- function(x = 10^seq(-6, 2, length.out = 100), ...) {

  y <- gpf_fraction_bound(x = x, ...)

  graphics::plot(x, y, log = "x", type = "l", ylab = expression(theta),
                 ylim = c(0, attr(y, "degree")))

  graphics::points(x[-1], diff(y) / max(diff(y)), type = "l", col = "gray")

  graphics::axis(side = 3, at = 1 / attr(y, "params")[1:attr(y, "degree")],
                 lwd = 0, lwd.ticks = 1)

}

#' Quick overview plots to demonstrate harmonic mean
#'
#' @noRd
gpf_fraction_micro_macro <- function(a, b, c = 0, x = 10^seq(-6, 2, length.out = 100),
                                     correlation = c(ab = 2, ac = 2, bc = 2, abc = 6)) {

  correlation <- c(correlation, c(ab = 0, bc = 0, ac = 0, abc = 0)[setdiff(
    c("ab", "bc", "ac", "abc"), names(correlation))])

  K1 = a + b + c
  ab = unname(correlation["ab"] * a * b)
  ac = unname(correlation["ac"] * a * c)
  bc = unname(correlation["bc"] * b * c)
  abc = unname(correlation["abc"] * a * b * c)
  K2 = (ab + ac + bc) / K1
  K3 = abc / (K2 * K1)

  macro <- c(K1 = K1, K2 = K2, K3 = K3)
  micro <- c(a = a, b = b, c = c, ab = ab, ac = ac, bc = bc, abc = abc)

  print(macro)
  print(micro)

  y_macro <- gpf_fraction_bound(x, binding_constants = macro, type = "macro")
  y_micro <- gpf_fraction_bound(x, binding_constants = micro, type = "micro")

  graphics::plot(x, y_macro, log = "x", type = "l",
                 ylab = expression(theta), ylim = c(0, 3 - (c == 0)),
                 col = "black")

  graphics::points(x, y_micro, type = "l", col = "red", lty = 2)

  for (i in micro[which(micro != 0)]) {

    graphics::points(x, sapply(x, gpf_fraction_bound, binding_constants = i, type = "micro"),
                     type = "l", lty = 2)

  }

}
#
# gpf_fraction_micro_macro(a = 1e-2, b = 1e2, c = 1e3)
# gpf_fraction_micro_macro(a = 1e3, b = 1e3, c = 0)

#' Fit binding isotherm
#'
#' @inheritParams fit_Kd
#' @inheritParams gpf_fraction_bound
#' @param INDEX A column variable that shall be used to split the data into at
#' most two sets, for which both will share higher order K. You will get other
#' estimates as .x and .y fits.
#' @param degree The order of the binding polynomial; if \code{NULL} it is inferred
#' from the number of variables on the LHS of \code{formula}. Up to third order
#' binding polynomials are supported.
#' @param correlation The correlation factors between the microscopic binding
#' constants; must be a named vector with "ab", "bc", "ac", "abc" as names.
#' @param start_K_d Upper and lower bounds of the K_d grid-start
#' parameters (see \link[nls.multstart:nls_multstart]{nls.multstart::nls_multstart}).
#'
#' @details
#' The \code{formula} LHS should evaluate to the observed binding isotherm, i.e.
#' the number of ligands bound at the given concentration (RHS). Column arithmetics
#' are accepted and even preferred as the number of variables in the LHS will be
#' used to determine the binding degree: \code{bound_1 + 2 * bound_2 ~ conc}
#' indicates a second degree binding polynom.
#'
#' If \code{type == "macro"}, the binding isotherm is computed for a macroscopic
#' (thermodynamic) process using \code{\link{gpf_macro}}.
#'
#' If \code{type == "micro"},
#' the intrinsic binding constants of the individual binding sites are computed
#' temptatively using correlations that yield to harmonic means. (In any system
#' this is non-sensical.) The vector must be strictly named with "ab", "bc", "ac",
#' "abc" as appropriate. Missing values will be set to 0, i.e., the result
#' will be the same as for \code{type == "macro"}.
#'
#' The grid start approach will take 540 points across all model estimates.
#'
#' If \code{INDEX} is a column name, the data is internally split into (exactly)
#' two groups, for which the higher order K are shared. Note that the grid start
#' approach will take 1296 points and computation may take longer.
#'
#' Up to third order binding polynoms can be computed.
#'
#' Note that all K_d values are internally log-transformed and returned as such
#' to allow for an equal search depth across all orders of magnitude using
#' \link[nls.multstart:nls_multstart]{nls.multstart::nls_multstart}.
#'
#' @examples
#' library(dplyr)
#' assay_file <- system.file("extdata", "gel_03.txt", package = "summerrband")
#' assay_data <- iqtl_meta(assay_file, list(conc = c(2^seq(10, 0), 0)))
#'
#' test3 <- assay_data %>%
#'   group_by(gel_id) %>%
#'   tidyr::pivot_wider(names_from = band_id, values_from = vol_frac,
#'   id_cols = c(conc, group_vars(.))) %>%
#'   model_cleanly_groupwise(fit_binding_isotherm, formula = band_1 + 2 * band_2 ~ conc,
#'                           newdata = data.frame(conc = 10^seq(-3, 3, length.out = 100)))
#'
#' tidyr::unnest(test3, tidy)
#'
#' library(ggplot2)
#' model_display(test3) + scale_x_log10()
#'
#' # fitting with shared estimates
#' assay_data_2 <- assay_data %>%
#'   group_by(gel_id) %>%
#'   tidyr::pivot_wider(names_from = band_id, values_from = vol_frac,
#'   id_cols = c(conc, group_vars(.))) %>%
#'   mutate(index_col = rep(c("a", "b"), times = 6))
#'
#' # these examples will take some computation time ...
#' #
#' # fit_binding_isotherm(assay_data_2, formula = band_1 + 2 * band_2 ~ conc, INDEX = index_col)
#' #
#' # test4 <- assay_data_2 %>%
#' #   model_cleanly_groupwise(fit_binding_isotherm, formula = band_1 + 2 * band_2 ~ conc,
#' #                           INDEX = index_col,
#' #                           newdata = data.frame(conc = 10^seq(-3, 3, length.out = 100)))
#'
#' @export
fit_binding_isotherm <- function(x, formula, degree = NULL, type = "macro", INDEX = NULL,
                                 limits_lower = c(-Inf, +Inf), limits_upper = c(-Inf, +Inf),
                                 limits_K_d = c(0, 1e3), start_K_d = 10^c(-1, 4),
                                 correlation = c(ab = 2, ac = 2, bc = 2, abc = 6)) {

  # generate observed binding isotherm according to formula virtually

  x <- dplyr::mutate(.data = x, .RL = !!formula.tools::lhs(formula))

  # determine degree of binding polynom to fit

  if (is.null(degree)) degree <- length(formula.tools::lhs.vars(formula))

  degree <- as.integer(degree)

  stopifnot(all(is.finite(degree), degree <= 3))  # currently only third degree implemented

  # construct formula to fit

  L0 <- formula.tools::rhs(formula)

  # proceed according to whether INDEX is a column variable or not

  INDEX <- rlang::enquo(INDEX)

  if (rlang::quo_is_null(INDEX)) {

    FML <- substitute(.RL ~ RL_isotherm(conc_L = L0, pK_d1, pK_d2, pK_d3, upper, lower,
                                        type = T0, correlation = C0),
                      list(L0 = L0, T0 = type, C0 = correlation))

    starts <- list(lower = c(0, 0.2), upper = c(0.5, 1.0),
                   pK_d1 = log10(start_K_d), pK_d2 = log10(start_K_d),
                   pK_d3 = log10(start_K_d))

    iters <- list(lower = 3, upper = 3, pK_d1 = 5, pK_d2 = 4, pK_d3 = 3)

    params <- list(lower = limits_lower, upper = limits_upper,
                   pK_d1 = log10(limits_K_d), pK_d2 = log10(limits_K_d),
                   pK_d3 = log10(limits_K_d))

  } else {

    INDEX <- rlang::as_name(INDEX)

    stopifnot(INDEX %in% colnames(x))

    FML <- substitute(.RL ~ RL_isotherm_shared(conc_L = L0, INDEX = I0,
                                               pK_d1.x, pK_d1.y, pK_d2, pK_d3,
                                               upper.x, lower.x, upper.y, lower.y,
                                               type = T0, correlation = C0),
                      list(L0 = L0, I0 = x[[INDEX]], T0 = type, C0 = correlation))

    starts <- list(lower.x = c(0, 0.2), upper.x = c(0.5, 1.0),
                   lower.y = c(0, 0.2), upper.y = c(0.5, 1.0),
                   pK_d1.x = log10(start_K_d), pK_d1.y = log10(start_K_d),
                   pK_d2 = log10(start_K_d), pK_d3 = log10(start_K_d))

    iters <- list(lower.x = 2, upper.x = 2,
                  lower.y = 2, upper.y = 2,
                  pK_d1.x = 3, pK_d1.y = 3,
                  pK_d2 = 3, pK_d3 = 3)

    params <- list(lower.x = limits_lower, upper.x = limits_upper,
                   lower.y = limits_lower, upper.y = limits_upper,
                   pK_d1.x = log10(limits_K_d), pK_d1.y = log10(limits_K_d),
                   pK_d2 = log10(limits_K_d), pK_d3 = log10(limits_K_d))

  }

  if (degree < 3) {

    iters$pK_d3 <- starts$pK_d3 <- params$pK_d3 <- NULL
    FML <- do.call("substitute", list(FML, list(pK_d3 = Inf)))

  }

  if (degree < 2) {

    iters$pK_d2 <- starts$pK_d2 <- params$pK_d2 <- NULL
    FML <- do.call("substitute", list(FML, list(pK_d2 = Inf)))

  }

  ll <- sapply(params, min, na.rm = TRUE, USE.NAMES = TRUE)
  ul <- sapply(params, max, na.rm = TRUE, USE.NAMES = TRUE)

  li <- sapply(starts, min, na.rm = TRUE, USE.NAMES = TRUE)
  ui <- sapply(starts, max, na.rm = TRUE, USE.NAMES = TRUE)

  eval(rlang::call2(.fn = "nls_multstart", .ns = "nls.multstart",
                    formula = stats::as.formula(FML),
                    data = x,
                    lower = ll, upper = ul,
                    iter = unlist(iters),
                    start_lower = li, start_upper = ui,
                    supp_errors = "Y"))

}

#' Function to fit an isotherm
#'
#' @noRd
RL_isotherm <- function(conc_L, pK_d1, pK_d2, pK_d3, upper, lower,
                        correlation = c(ab = 1, bc = 1, ac = 1, abc = 1),
                        type = "macro") {

  # use log-transformed K_d to allow equal search space of nls_multstart

  if (type == "micro") {

    correlation <- c(correlation, c(ab = 0, bc = 0, ac = 0, abc = 0)[setdiff(
      c("ab", "bc", "ac", "abc"), names(correlation))])

    a <- 10^(-pK_d1)
    b <- 10^(-pK_d2)
    c <- 10^(-pK_d3)
    ab <- unname(correlation["ab"] * a * b)
    ac <- unname(correlation["ac"] * a * c)
    bc <- unname(correlation["bc"] * b * c)
    abc <- unname(correlation["abc"] * a * b * c)

    params <- c(a = a, b = b, c = c, ab = ab, ac = ac, bc = bc, abc = abc)

  } else {

    params <- c(K1 = 10^(-pK_d1), K2 = 10^(-pK_d2), K3 = 10^(-pK_d3))

  }

  lower + (upper - lower) * gpf_fraction_bound(
    x = conc_L, binding_constants = params, type = type)

}

#' Function to fit an isotherm with shared pK_d2 and pK_d3
#'
#' @noRd
RL_isotherm_shared <- function(conc_L, INDEX = NULL,
                               pK_d1.x, pK_d1.y = pK_d1.x, pK_d2, pK_d3,
                               upper.x, lower.x,
                               upper.y = upper.x, lower.y = lower.y,
                               type = "macro", correlation = c(ab = 0, bc = 0, ac = 0, abc = 0)) {

  INDEX <- as.factor(INDEX)

  if (length(levels(INDEX)) == 2) {

    conc_L <- split(conc_L, INDEX, drop = FALSE)

    rlx <- RL_isotherm(conc_L[[1]], pK_d1.x, pK_d2, pK_d3, upper.x, lower.x,
                       type = type, correlation = correlation)
    rly <- RL_isotherm(conc_L[[2]], pK_d1.y, pK_d2, pK_d3, upper.y, lower.y,
                       type = type, correlation = correlation)

    res <- list(rlx, rly)
    names(res) <- names(conc_L)

    res <- unsplit(res, INDEX, drop = TRUE)

  } else {

    if(length(levels(INDEX)) > 2) warning("Data splits into more than 2 groups;
                                             not grouping at all now.")

    res <- RL_isotherm(conc_L, pK_d1.x, pK_d2, pK_d3, upper.x, lower.x,
                       type = type, correlation = correlation)

  }

  res

}

#' Augmenting function for shared isotherm fits
#'
#' @param x A fitted model returned from \code{\link{fit_binding_isotherm}}
#' @param newdata The new data to fit; must use the same column names that
#' were used to create the model.
#' @inheritParams fit_binding_isotherm
#'
#' @export
augment_shared_isotherms <- function(x, newdata, INDEX = NULL) {

  INDEX <- rlang::enquo(INDEX)
  INDEX <- rlang::as_name(INDEX)

  stopifnot(INDEX %in% colnames(newdata))

  new.call <- as.list(formula.tools::rhs(x$m$formula()))
  new.call$INDEX <- newdata[[INDEX]]
  new.call <- as.call(new.call)

  cbind(newdata, .fitted = with(newdata, eval(new.call, envir = as.list(coef(x)))))

}

