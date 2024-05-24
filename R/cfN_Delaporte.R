#' @title Characteristic function of Delaporte distribution
#'
#' @description
#' \code{cfN_Delaporte(t, a, b, c)} evaluates the characteristic function \eqn{cf(t)} of the
#' Delaporte distribution with the parameters \code{a} (parameter of variable mean, \code{a > 0})
#' \code{b} (parameter of variable mean, \code{b > 0} ), and \code{c} (fixed mean, \code{c > 0}), i.e.
#' \deqn{cfN_Delaporte(t, a, b, c) = (b/(1+b))^a * (1-e^(1i*t)/(b+1))^(-a) * exp(-c*(1-e^(1i*t))).}
#' For more details see [4].
#'
#' \code{cfN_Delaporte(t, a, b, c, cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = cfN_Delaporte(-1i*log(cfX(t)), a, b, c),} where \code{cfX} is function
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
#' \url{https://en.wikipedia.org/wiki/Delaporte_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param a variable mean, \code{a > 0}.
#' @param b variable mean, \code{b > 0}.
#' @param c fixed mean, \code{c > 0}.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Delaporte distribution.
#'
#' @note Ver.: 16-Sep-2018 18:58:36 (consistent with Matlab CharFunTool v1.3.0, 15-Nov-2016 13:36:26).
#'
#' @example R/Examples/example_cfN_Delaporte.R
#'
#' @export
#'
cfN_Delaporte <- function(t, a, b, c, cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t) || missing(a) || missing(b) || missing(c)) {
    stop("Enter input parameters t, a, b, c.")
  }

  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {
    expit <- exp(1i * t)
  } else {
    expit = cfX(t)
  }

  cf <- (b / (1 + b)) ^ a * (1 - expit / (b + 1)) ^ (-a) * exp(-c * (1 -
                                                                       expit))

  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
