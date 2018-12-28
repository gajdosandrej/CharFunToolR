#' @title Characteristic function of a GAUSSIAN (standard normal) distribution
#'
#' @description
#' \code{cfS_Gaussian(t, mu, sigma, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of a GAUSSIAN (standard normal) distribution.
#'
#' cfS_Gaussian is an ALIAS of the more general function cf_Normal, used
#' to evaluate the characteristic function of a linear combination
#' of independent normally distributed random variables.
#'
#' The characteristic function of the standard normally distributed random variable, \eqn{X ~ N(0,1)},
#' is defined by \deqn{cf(t) = cfS_Gaussian(t) = exp(-t^2/2).}
#'
#' @family Continuous Probability Distribution
#' @family Symmetric Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param mu vector of mean values parameters.
#' @param sigma vector of dispersion parameters.
#' @param coef vector of coefficients of the linear combination of Normally distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient.
#'
#' @return Characteristic function \eqn{cf(t)} of a GAUSSIAN (standard normal) distribution.
#'
#' @note Ver.: 16-Sep-2018 19:08:06 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cfS_Gaussian.R
#'
#' @export
#'
cfS_Gaussian <- function(t, mu, sigma, coef, niid) {
  if (missing(mu)) {
    mu <- vector()
  }

  if (missing(sigma)) {
    sigma <- vector()
  }

  if (missing(coef)) {
    coef <- vector()
  }

  if (missing(niid)) {
    niid <- vector()
  }

  cf <- cf_Normal(t, mu, sigma, coef, niid)

  return(cf)

}
