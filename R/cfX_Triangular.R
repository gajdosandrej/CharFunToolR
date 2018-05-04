#' @title Characteristic function of Triangular distribution
#'
#' @description
#' \code{cfX_Triangular(t, a, b, c)} (\code{a} \eqn{\le} \code{c} \eqn{\le} \code{b})
#' evaluates the characteristic function \eqn{cf(t)} of the Triangular distribution
#' on the interval \eqn{(}\code{a}, \code{b}\eqn{)} with mode \code{c}
#' (Triangular distribution with \eqn{mean = a + b + c)/3} and \eqn{variance = 1/18(a^2 + b^2 + c^2 + - ab - ac - bc)}).
#' \deqn{cfX_Triangula(t, a, b, c) = -2((b-c)exp(iat) - (b-a)exp(ict) + (c-a)exp(ibt))/((b-a)(c-a)(b-c)t^2).}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Triangular_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param a number, default value \eqn{a = -1}.
#' @param b number, default value \eqn{b = 1}.
#' @param c number, (\eqn{a \le c \le b}), default value \eqn{c = 0}
#' @return Characteristic function \eqn{cf(t)} of the Triangular distribution.
#'
#' @example R/Examples/example_cfX_Triangular.R
#'
#' @export
#'
cfX_Triangular <- function(t,
                           a = -1,
                           b = 1,
                           c = 0) {
  szt <- dim(t)
  t <- c(t)

  cf <-
    -2 * ((b - c) * exp(1i * a * t) - (b - a) * exp(1i * c * t) + (c - a) *
            exp(1i * b * t)) / ((b - a) * (c - a) * (b - c) * t ^ 2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
