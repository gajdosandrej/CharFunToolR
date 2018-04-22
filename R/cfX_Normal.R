#' @title Characteristic function of Normal distribution
#'
#' @description
#' \code{cfX_Normal(t, mean, variance)} evaluates the characteristic function \eqn{cf(t)} of
#' the Normal distribution with \eqn{mean = mean} and \eqn{variance = variance}: \eqn{N(mean, variance)};\cr
#' \eqn{cfX_Normal(t, mean, variance) = exp(imeant -1/2variance^2t^2)}.
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}.
#'
#' @param t numerical values (number, vector...)
#' @param mean real number, mean or expextation of the distribution, default value \eqn{mean = 0}
#' @param variance real number, standard deviation, \eqn{variance > 0}, default value \eqn{variance = 1}
#' @return Characteristic function \eqn{cf(t)} of the normal distribution.
#'
#' @example R/Examples/example_cfX_Normal.R
#'
#' @export
#'
cfX_Normal <- function(t, mean = 0, variance = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- exp(1i * mean * t - 1 / 2 * variance ^ 2 * t ^ 2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
