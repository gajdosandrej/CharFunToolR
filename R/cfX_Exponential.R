#' @title Characteristic function of the EXPONENTIAL distribution
#'
#' @description
#' \code{cfX_Exponential(t, lambda)} evaluates characteristic function of the EXPONENTIAL distribution
#' with the rate parameter \code{lambda > 0}.
#'
#' \code{cfX_Exponential} is an ALIAS NAME of the more general function \code{cf_Exponential},
#' used to evaluate the characteristic function of a linear combination
#'
#' of independent EXPONENTIAL distributed random variables.
#' The characteristic function of the EXPONENTIAL distribution is
#' \deqn{cf(t) = \lambda / (\lambda - 1i*t).}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Exponential_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param lambda vector of the 'rate' parameters \code{lambda > 0}. If empty, default value is \code{lambda = 1}.
#' @param coef vector of coefficients of the linear combination of Exponentially distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient.
#'
#' @return Characteristic function \eqn{cf(t)} of the EXPONENTIAL distribution.
#'
#' @note Ver.: 16-Sep-2018 19:21:37 (consistent with Matlab CharFunTool v1.3.0, 24-Jun-2017 10:07:43).
#'
#' @example R/Examples/example_cfX_Exponential.R
#'
#' @export
#'
cfX_Exponential <- function(t, lambda, coef, niid) {
  if (missing(lambda)) {
    lambda <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  cf <- cf_Exponential(t, lambda, coef, niid)

  return(cf)
}
