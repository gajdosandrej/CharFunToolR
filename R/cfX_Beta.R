#' @title Characteristic function of the BETA distribution
#'
#' @description
#' \code{cfX_Beta(t, alpha, beta)} evaluates the characteristic function \eqn{cf(t)} of
#' the Beta distribution with the parameter \code{alpha} (shape, \eqn{\alpha > 0}) and \code{beta} (shape, \eqn{\beta > 0})
#' defined on the interval \eqn{(0,1)}, i.e. beta distribution with the \eqn{Mean = \alpha / (\alpha + \beta)}
#' and the \eqn{Variance = (\alpha*\beta) / ((\alpha+\beta)^2*(\alpha+\beta+1))}.
#' Then, the standard deviation is given by \eqn{STD = sqrt(Variance)}
#' i.e.
#' \deqn{cf(t) = cfX_Beta(t,\alpha,\beta) = 1F1(\alpha ,\alpha + \beta , i*t),}
#' where \eqn{1F1(.;.;.)} is the Confluent hypergeometric function.
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha shape, \code{alpha > 0}, default value \code{alpha = 1}.
#' @param beta shape, \code{beta > 0}, default value \code{beta = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of the BETA distribution.
#'
#' @example R/Examples/example_cfX_Beta.R
#'
#' @export
#'
cfX_Beta <- function(t, alpha = 1, beta = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- hypergeom1F1(1i * t, alpha, alpha + beta)


  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
