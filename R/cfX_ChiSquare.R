#' @title Characteristic function of the CHI-SQUARE distribution
#'
#' @description
#' \code{cfX_ChiSquare(t, df, ncp, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of the CHI-SQUARE distribution with \code{df > 0} degrees of freedom and the non-cetrality parameter \code{ncp > 0}.
#'
#' \code{cfX_ChiSquare} is an ALIAS NAME of the more general function
#' \code{cf_ChiSquare}, used to evaluate the characteristic function of a linear combination
#' of independent CHI-SQUARE distributed random variables.
#'
#' The characteristic function of the CHI-SQUARE distribution is defined
#' by \deqn{cf(t) = (1 - 2*i*t )^(df/2) * exp((i*t*ncp)/(1-2*i*t)).}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Chi-squared_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param df the degrees of freedom parameter \code{df > 0}. If empty, default value is \code{df = 1}.
#' @param npc the non-centrality parameter \code{ncp > 0}. If empty, default value is \code{ncp = 0}.
#'
#' @return Characteristic function \eqn{cf(t)} of the CHI-SUQARE distribution.
#'
#' @example R/Examples/example_cfX_ChiSquare.R
#'
#' @export
#'
cfX_ChiSquare <- function(t, df, ncp, coef, niid) {
  if (missing(df)) {
    df <- vector()
  }
  if (missing(ncp)) {
    ncp <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- vector()
  }

  cf <- cf_ChiSquare(t, df, ncp, coef, niid)

  return(cf)
}
