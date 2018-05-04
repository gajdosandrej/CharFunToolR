#' @title Characteristic function of the GAMMA distribution
#'
#' @description
#' \code{cfX_Gamma(t, alpha, beta, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of the GAMMA distribution with the shape parameter \code{alpha > 0} and the rate parameter \code{beta > 0}.
#'
#' \code{cfX_Gamma} is an ALIAS NAME of the more general function \code{cf_Gamma},
#' used to evaluate the characteristic function of a linear combination
#' of independent GAMMA distributed random variables.
#'
#' The characteristic function of the GAMMA distribution is defined by
#' \deqn{cf(t) = (1 - i*t/\beta)^(-\alpha).}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Gamma_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha the shape parameter \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta the rate (\eqn{1/scale}) parameter \code{beta > 0}. If empty, default value is \code{beta = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of the GAMMA distribution.
#'
#' @example R/Examples/example_cfX_Gamma.R
#'
#' @export
#'
cfX_Gamma <- function(t, alpha, beta, coef, niid) {
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
    niid <- vector()
  }

  cf <- cf_Gamma(t, alpha, beta, coef, niid)

  return(cf)
}
