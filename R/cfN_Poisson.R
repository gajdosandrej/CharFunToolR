#' @title Characteristic function of the Poisson distribution
#'
#' @description
#' \code{cfN_Poisson(t, lambda, cfX)} evaluates the characteristic function \eqn{cf(t)} of the
#' Poisson distribution with the rate parameter \code{lambda > 0}, i.e.
#' \deqn{cfN_Poisson(t, \lambda) = exp(\lambda*(exp(1i*t)-1))}.
#' For more details see [4].
#'
#' \code{cfN_Poisson(t, lambda, cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = cfN_Poisson(-1i*log(cfX(t)), \lambda),} where \code{cfX} is function
#' handle of the characteristic function \eqn{cfX(t)} of a continuous distribution
#' and/or random variable \eqn{X}.
#'
#' Note that such CF is characteristic function of the compound distribution,
#' i.e. distribution of the random variable \eqn{Y = X_1 + ... + X_N}, where \eqn{X_i ~ F_X}
#' are i.i.d. random variables with common CF \eqn{cfX(t)}, and \eqn{N ~ F_N} is
#' independent RV with its CF given by \eqn{cfN(t)}.
#'
#' @family Discrete Probability Distribution
#'
#' @references
#' [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
#'     distribution based on numerical inversion of the compound empirical
#'     characteristic function of frequency and severity. Preprint submitted
#'     to Insurance: Mathematics and Economics.
#'
#' [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM
#'     for computing distributions of collective risk models. Preprint submitted
#'     to Journal of Statistical Software.
#'
#' [3] WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#'     An alternative tool to form the probability distribution
#'     of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
#'     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN 3-87773-025-6.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Poisson_distribution}, \cr
#' \url{https://en.wikipedia.org/wiki/Compound_Poisson_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param lambda rate, \code{lambda > 0}, default value \code{lambda = 1}.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Poisson distribution.
#'
#' @example R/Examples/example_cfN_Poisson.R
#'
#' @export
#'
cfN_Poisson <- function(t, lambda = 1, cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t)) {
    stop("Enter input parameter t.")
  }

  ## Characteristic function of the (compound) Poisson distribution
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {
    expit <- exp(1i * t)
  } else {
    expit = cfX(t)
  }

  cf <- exp(lambda * (expit - 1))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
