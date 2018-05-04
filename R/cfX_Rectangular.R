#' @title Characteristic function of Rectangular distribution
#'
#' @description
#' \code{cfX_Rectangular(t, a, b)} evaluates the characteristic function \eqn{cf(t)} of
#' the Rectangular distribution on the interval \eqn{(}\code{a}, \code{b}\eqn{)}
#' (Rectangular distribution with \eqn{mean = (a + b)/2} and \eqn{variance = 1/12(b - a)^2)}
#' \eqn{cfX_Rectangular(t, a, b) = (exp(ibt) - exp(iat))/(i(b - a)t)}.
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param a number, default value \eqn{a = -1}.
#' @param b number, default value \eqn{b = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of the Rectangular distribution.
#'
#' @example R/Examples/example_cfX_Rectangular.R
#'
#' @export
#'
cfX_Rectangular <- function(t, a = -1, b = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- (exp((0 + 1i) * b * t) - exp((0 + 1i) * a * t)) / ((0 + 1i) * (b - a) *
                                                             t)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
