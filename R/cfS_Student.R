#' @title Characteristic function of the STUDENT's t-distribution
#'
#' @description
#' \code{cfS_Student(t, df, mu, sigma, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of the STUDENT's t-distribution with \code{df > 0} degrees of freedom.
#'
#' \code{cfS_Student} is an ALIAS of the more general function \code{cf_Student}, used
#' to evaluate the characteristic function of a linear combination
#' of independent (location-scale) STUDENT's t-distributed random variables.
#'
#' The characteristic function of the STUDENT's t-distribution with \eqn{df} degrees of freedom is defined by
#' \deqn{cf(t) = cfS_Student(t,df) = besselk(df/2,abs(t)*sqrt(df),1) * exp(-abs(t)*sqrt(df)) * (sqrt(df)*abs(t))^(df/2) / 2^(df/2-1)/gamma(df/2).}
#'
#' @family Continuous Probability distribution
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic
#' function: An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Student's_t-distribution}.
#'
#' @param t  vector or array of real values, where the CF is evaluated.
#' @param df the degrees of freedom, \code{df > 0}. If empty, the default value is \code{df = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of the STUDENT's t-distribution.
#'
#' @example R/Examples/example_cfS_Student.R
#'
#' @export
#'
cfS_Student <- function(t, df = 1, mu, sigma, coef, niid) {
  if (missing(df)) {
    df <- vector()
  }
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
    niid <- numeric()
  }

  cf <- cf_Student(t, df, mu, sigma, coef, niid)

  return(cf)

}
