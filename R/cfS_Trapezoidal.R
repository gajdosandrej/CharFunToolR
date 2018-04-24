#' @title Characteristic function of the zero-mean symmetric TRAPEZOIDAL distribution
#'
#' @description
#' \code{cfS_Trapezoidal(t, lambda, coef, niid)} evaluates
#' the characteristic function of the zero-mean symmetric TRAPEZOIDAL distribution defined on the interval \eqn{(-1,1)}.
#'
#' \code{cfS_Trapezoidal} is an ALIAS of the more  general function
#' \code{cf_TrapezoidalSymmetric}, used to evaluate the characteristic function
#' of a linear combination of independent TRAPEZOIDAL distributed random variables.
#'
#' The characteristic function of \eqn{X ~ TrapezoidalSymmetric(\lambda)}, where
#' \eqn{0\le \lambda \le 1} is the offset parameter is defined by
#' \deqn{cf(t) = (sin(w*t)/(w*t))*(sin((1-w)*t)/((1-w)*t))},
#' where \eqn{w  = (1+\lambda)/2}.
#'
#' @family Continuous Probability distribution
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Trapezoidal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param lambda parameter of the offset, \eqn{0 \le} \code{lambda} \eqn{\le 1}.
#' If empty, default value is \code{lambda = 0}.
#'
#' @return Characteristic function \eqn{cf(t)} of the zero-mean symmetric TRAPEZOIDAL distribution.
#'
#' @example R/Examples/example_cfS_Trapezoidal.R
#'
#' @export
#'
cfS_Trapezoidal <- function(t, lambda, coef, niid) {
  if (missing(lambda)) {
    lambda <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  cf <- cf_TrapezoidalSymmetric(t, lambda, coef, niid)

  return(cf)
}

# cfS_Trapezoidal <- function(t, a = 1, c = 1/3) {
#   szt <- dim(t)
#   t <- c(t)
#
#   w = (1+c/a)/2
#
#   cf <- cfS_Rectangular(w*a*t)*cfS_Rectangular((1-w*a)*t)
#     cf[t == 0] <- 1
#
#   dim(cf) <- szt
#
#   return(cf)
# }
