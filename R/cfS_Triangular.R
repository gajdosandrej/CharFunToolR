#' @title Characteristic function of the zero-mean symmetric TRIANGULAR distribution
#'
#' @description
#' \code{cfS_Triangular(t, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of the zero-mean symmetric TRIANGULAR distribution defined on the interval \eqn{(-1,1)}.
#'
#' \code{cfS_Triangular} is an ALIAS of the more general function \code{cf_TriangularSymmetric},
#' used to evaluate the characteristic function of a linear combination
#' of independent TRIANGULAR distributed random variables.
#'
#' The characteristic function of \eqn{X ~ TriangularSymmetric} is defined by
#' \deqn{cf(t) = (2-2*cos(t))/t^2.}
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Triangular_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#'
#' @return Characteristic function \eqn{cf(t)} of the zero-mean symmetric TRIANGULAR distribution.
#'
#' @example R/Examples/example_cfS_Triangular.R
#'
#' @export
#'
cfS_Triangular <- function(t, coef, niid) {
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  cf <- cf_TriangularSymmetric(t, coef, niid)

  return(cf)
}
