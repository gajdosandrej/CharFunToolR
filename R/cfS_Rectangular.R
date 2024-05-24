#' @title Characteristic function of the zero-mean symmetric RECTANGULAR distribution
#'
#' @description
#' \code{cfS_Rectangular(t, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of the zero-mean symmetric RECTANGULAR distribution defined on the interval \eqn{(-1,1)}.
#'
#' \code{cfS_Rectangular} is an ALIAS of the more general function
#' \code{cf_RectangularSymmetric}, used to evaluate the characteristic function
#' of a linear combination of independent RECTANGULAR distributed random variables.
#'
#' The characteristic function of \eqn{X ~ RectangularSymmetric} is defined by
#' \deqn{cf(t) = sinc(t) = sin(t)/t.}
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
#' \url{https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param coef vector of coefficients of the linear combination of Rectangular distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient.
#'
#' @return Characteristic function \eqn{cf(t)} of the zero-mean symmetric RECTANGULAR distribution.
#'
#' @note Ver.: 16-Sep-2018 19:09:05 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cfS_Rectangular.R
#'
#' @export
#'
cfS_Rectangular <- function(t, coef, niid) {
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  cf <- cf_RectangularSymmetric(t, coef, niid)

  return(cf)

}
