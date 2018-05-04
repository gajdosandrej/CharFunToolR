#' @title Characteristic function of the INVERSE GAMMA distribution
#'
#' @description
#' \code{cfX_InverseGamma(t, alpha, beta, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of the INVERSE GAMMA distribution with the shape parameter \code{alpha > 0} and the rate parameter \code{beta > 0}.
#'
#' \code{cfX_InverseGamma} is an ALIAS NAME of the more general function \code{cf_InverseGamma},
#' used to evaluate the characteristic function of a linear combination
#' of independent INVERSE GAMMA distributed random variables.
#'
#' The characteristic function of the GAMMA distribution is defined
#' by \deqn{cf(t) =  2 / gamma(\alpha) * (-1i*\beta*t).^(\alpha/2) * besselk(\alpha,sqrt(-4i*\beta*t)).}
#'
#' @family Continuous Probability Distribution
#'
#' @references
#' WITKOVSKY, V.: Computing the distribution of a linear combination
#' of inverted gamma variables, Kybernetika 37 (2001), 79-90.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha the shape parameter \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta the rate (1/scale) parameter \code{beta > 0}. If empty, default value is \code{beta = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of the INVERSE GAMMA distribution.
#'
#' @example R/Examples/example_cfX_InverseGamma.R
#'
#' @export
#'
cfX_InverseGamma <- function(t, alpha, beta, coef, niid) {
  if (missing(alpha)) {
    alpha <- vector()
  }
  if (missing(beta)) {
    beta <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  cf <- cf_InverseGamma(t, alpha, beta, coef, niid)

  return(cf)
}
