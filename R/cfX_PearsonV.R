#' @title Characteristic function of Pearson type V distribution
#'
#' @description
#' \code{cfX_PearsonV(t, alpha, beta)} evaluates the characteristic function
#' of the Pearson type V distribution with the parameters \eqn{alpha} (shape, \code{alpha > 0}) and
#' \eqn{beta} (scale, \code{beta > 0}), computed for real vector argument \code{t}, i.e.
#' \deqn{cfX_PearsonV(t, \alpha, \beta) = (2/gamma(\alpha)) * (-1i*t/\beta)^(\alpha/2) * besselk(\alpha,2*sqrt(-1i*t/\beta)),}
#' where \eqn{besselk(a,z)} denotes the modified Bessel function of the second order.
#'
#' @family Continuous Probability distribution
#'
#'  @references
#' [1] WITKOVSKY, V.: On the exact computation of the density and
#' of the quantiles of linear combinations of t and F random variables.
#' Journal of Statistical Planning and Inference 94 (2001), 1-13.
#'
#' [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution of output quantity
#' in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' [3] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss distribution
#' based on numerical inversion of the compound empirical characteristic function
#' of frequency and severity. Working Paper. Insurance: Mathematics and Economics.
#'
#' [4] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for computing distributions
#' of collective risk models.  Working Paper. Journal of Statistical Software.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Pearson_distribution}.
#'

#' @param t vector of real values where the CF is evaluated, i.e. \eqn{CF(t)}.
#' @param alpha scalar shape parameter of the Pearson VI distribution, \code{alpha > 0}.
#' @param beta scalar scale parameter of the Pearson VI distribution, \code{beta > 0}.
#'
#' @return Characteristic function \eqn{cf(t)} of the Pearson type V distribution.
#'
#' @example R/Examples/example_cfX_PearsonV.R
#'
#' @export
#'
cfX_PearsonV <- function(t, alpha = 1, beta = 1) {
  szt <- dim(t)

  cf <-
    unlist(lapply(t, function(t)
      tryCatch(
              Bessel::BesselK(2 * sqrt((0 - 1i * t) / beta), nu = alpha),
        error = function(e)
          0
      )))
  cf <- (2 / gamma(alpha)) * ((0 - 1i) * t / beta) ^ (alpha / 2) * cf

  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
