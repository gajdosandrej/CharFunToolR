#' @title Characteristic function of the symmetric zero-mean Arcsine distribution
#'
#' @description
#' \code{cfS_Arcsine(t, coef, niid)} evaluates the characteristic function \eqn{cf(t)} of
#' the zero-mean symmetric Arcsine distribution on the interval
#' \eqn{(-1,1)}.
#'
#' \code{cfS_Arcsine} is an ALIAS of the more general function
#' \code{cf_ArcsineSymmetric}, used to evaluate the characteristic function of a
#' linear combination of independent ARCSINE distributed random variables.
#'
#' The characteristic function of the symmetric ARCSINE distribution is \eqn{cf(t) = besselj(0,t)}.
#'
#' @family Continuous Probability Distribution
#' @family Symmetric Probability Distribution
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic
#' function: An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Arcsine_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param coef vector of coefficients of the linear combination of Arcsine distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient.
#'
#' @return Characteristic function \eqn{cf(t)} of the Arcsine distribution.
#'
#' @note Ver.: 16-Sep-2018 19:06:57 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cfS_Arcsine.R
#'
#' @export
#'
cfS_Arcsine <- function(t, coef, niid) {
  if (missing(coef)) {
    coef <- vector()
  }

  if (missing(niid)) {
    niid <- vector()
  }

  cf <- cf_ArcsineSymmetric(t, coef, niid)

  return(cf)
}
