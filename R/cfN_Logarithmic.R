#' @title Characteristic function of the Logarithmic distribution
#'
#' @description
#' \code{cfN_Logarithmic(t, p, cfX)} evaluates the characteristic function \eqn{cf(t)} of the
#' Logarithmic distribution defined on non-negative integers \eqn{n =0, 1, \ldots},
#' with the parameters \code{p} (success probability \code{p} in \eqn{[0,1]}), i.e.
#' \deqn{cfN_Logarithmic(t, p) = log(1 - p * exp(1i*t)) / log(1 - p).}
#' For more details see [4].
#'
#' \code{cfN_Logarithmic(t,p,cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = cfN_Logarithmic(-1i*log(cfX(t)), p),} where \code{cfX} is function
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
#' \url{https://en.wikipedia.org/wiki/Logarithmic_distribution}.
#'
#' @param t numerical values (number, vector...).
#' @param p success probability, \eqn{0 \le} \code{p} \eqn{\le 1}, default value \code{p = 1}.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Logarithmic distribution.
#'
#' @note Ver.: 16-Sep-2018 19:01:18 (consistent with Matlab CharFunTool v1.3.0, 15-Nov-2016 13:36:26).
#'
#' @example R/Examples/example_cfN_Logarithmic.R
#'
#' @export
#'
cfN_Logarithmic <- function(t, p = 1, cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t)) {
    stop("Enter input parameter t.")
  }

  ## Characteristic function of the (compound) Logarithmic distribution
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {
    expit <- exp(1i * t)
  } else {
    expit = cfX(t)
  }

  cf <- log(1 - p * expit) / log(1 - p)

  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
