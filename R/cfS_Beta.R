#' @title Characteristic function of the zero-mean symmetric BETA distribution
#'
#' @description
#' \code{cfS_Beta(t, theta, coef, niid)} evaluates the characteristic function \eqn{cf(t)} of
#' the zero-mean symmetric BETA distribution defined on the interval \eqn{(-1,1)}.
#'
#' \code{cfS_Beta} is an ALIAS of the more general function \code{cf_BetaSymmetric},
#' used to evaluate the characteristic function of a linear combination
#' of independent BETA distributed random variables.
#'
#' The characteristic function of \eqn{X ~ BetaSymmetric(\theta)} is defined by
#' \deqn{cf(t) = cf_BetaSymmetric(t,\theta) = gamma(1/2+\theta) * (t/2)^(1/2-\theta) * besselj(\theta-1/2,t).}
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
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param theta the 'shape' parameter \code{theta > 0}. If empty, default value is \code{theta = 1}.
#' @param coef vector of coefficients of the linear combination of Beta distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient.
#'
#' @return Characteristic function \eqn{cf(t)} of the Beta distribution.
#'
#' @note Ver.: 16-Sep-2018 19:07:26 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cfS_Beta.R
#'
#' @export
#'
cfS_Beta <- function(t, theta = 1, coef, niid) {
  if (missing(theta)) {
    theta <- vector()
  }

  if (missing(coef)) {
    coef <- vector()
  }

  if (missing(niid)) {
    niid <- numeric()
  }

  cf <- cf_BetaSymmetric(t, theta, coef, niid)

  return(cf)

}
