#' @title Characteristic function of the Polya-Eggenberger distribution
#'
#' @description
#' \code{cfN_PolyaEggenberger(t, a, b, m)} evaluates the characteristic function \eqn{cf(t)} of the
#' Polya-Eggenberger distribution with the parameters \code{a} (\code{a} real), \code{b} (\code{b} real),
#' and \code{m} (\code{m} integer), i.e.
#' \deqn{cfN_PolyaEggenberger(t, a, b, m) = 2F1(-m,a,a+b,1-e^(1i*t)),}
#' where \eqn{2F1} denotes the Gauss hypergeometric function. For more details see [4], p. 525.
#'
#' \code{cfN_PolyaEggenberger(t, a, b, m, cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = cfN_PolyaEggenberge(-1i*log(cfX(t)), a, b, m),} where \code{cfX} is function
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
#' @param t vector or array of real values, where the CF is evaluated.
#' @param a real number.
#' @param b real number.
#' @param m integer.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Polya-Eggenberger distribution.
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
#' @example R/Examples/example_cfN_PolyaEggenberger.R
#'
#' @export
#'
cfN_PolyaEggenberger <- function(t, a, b, m, cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t) || missing(a) || missing(b) || missing(m)) {
    stop("Enter input parameters t, a, b, m.")
  }

  ## Characteristic function of the (compound) Polya-Eggenberger distribution
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {
    expit <- exp(1i * t)
  } else {
    expit = cfX(t)
  }

  const1 <- 1
  const2 <- rep.int(1, length(t))
  cf <- const2

  for (i in 0:(m - 1)) {
    const1 <- const1 * (b + i) / (a + b + i)
    const2 <-
      (-m + i) * (a + i) / (-m - b + 1 + i) / (i + 1) * const2 * expit
    cf <- cf + const2
  }

  cf = cf * const1
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
