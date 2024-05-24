#' @title Characteristic function of Generalized Poisson distribution
#'
#' @description
#' \code{cfN_GeneralizedPoisson(t, a, p, cfX)} evaluates the characteristic function \eqn{cf(t)} of the
#' Generalized-Poisson distribution, with the parameters \code{a} (variable mean, \eqn{a > 0}),
#' and \code{p} (success probability, \eqn{0 \le p \le 1}) i.e.
#' \deqn{cfN_GeneralizedPoisson(t, a, p) = exp(a*(sum_{j=1}^Inf ((p*j)^(j-1)*e^(-p*j)/j!)*e^(1i*t*j)-1)).}
#'For more details see [4], p. 93.
#'
#' The Generalized-Poisson distribution is equivalent
#' with the Borel-Tanner distribution with parameters (p,m),
#' see WIMMER & ALTMANN (1999), where \eqn{m ~ Poisson(a)}.
#'
#' \code{cfN_GeneralizedPoisson(t, a, p, cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = cfN_GeneralizedPoisson(-1i*log(cfX(t)), a, p),} where \code{cfX} is function
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
#' @param t vector or array of real values, where the CF is evaluated.
#' @param a variable mean, \code{a > 0}.
#' @param p success probability, \eqn{0 \le} \code{p} \eqn{\le 1}, default value \code{p = 1/2}.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Poisson distribution.
#'
#' @note Ver.: 16-Sep-2018 18:59:23 (consistent with Matlab CharFunTool v1.3.0, 15-Nov-2016 13:36:26).
#'
#' @example R/Examples/example_cfN_GeneralizedPoisson.R
#'
#' @export
#'
cfN_GeneralizedPoisson <- function(t, a, p, cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t) || missing(a) || missing(p)) {
    stop("Enter input parameters t, a, p.")
  }

  ## Characteristic function of the Generalized-Poisson distribution
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {
    expit <- exp(1i * t)
  } else {
    expit = cfX(t)
  }

  exppt = exp(-p) * expit
  c     = exppt

  cf    = c

  N     = ceiling(-36.84 / (1 + log(p) - p))
  for (j in seq(1:N))
  {
    c  = p * ((1 + j) / j) ^ (j - 1) * exppt * c
    cf = cf + c
  }

  cf = exp(a * (cf - 1))

  cf[t == 0] = 1


  dim(cf) <- szt

  return(cf)
}
