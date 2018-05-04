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
#'
#' @return Characteristic function \eqn{cf(t)} of the EXPONENTIAL distribution.
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
