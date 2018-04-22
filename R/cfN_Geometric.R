#' @title Characteristic function of the Geometric distribution
#'
#' @description
#' \code{cfN_Geometric(t, p, type, cfX)} evaluates the characteristic function \eqn{cf(t)} of the
#' Geometric distribution.
#'
#'  The standard Geometric distribution (type = "standard" or "zero") is
#'  defined on non-negative integers \eqn{k = 0,1, \ldots} .
#'
#'  The shifted Geometric distribution (type = "shifted")
#'  is defined on positive integers \eqn{k = 1,2, \ldots} .
#'
#' Both types are parametrized by the success probability parameter \code{p} in \eqn{[0,1]}), i.e.
#' \eqn{cfN_Geometric(t, p, "standard") = p / (1 - (1-p) * exp(1i*t))},
#' \eqn{cfN_Geometric(t, p, "shifted")  = exp(1i*t) * (p / (1 - (1-p) * exp(1i*t)))}.
#' For more details see [4].
#'
#' \code{cfN_Geometric(t, p, type, cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = Geometric(-1i*log(cfX(t)), p),}
#' where \code{cfX} is function
#' handle of the characteristic function cfX(t) of a continuous distribution
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
#' \url{https://en.wikipedia.org/wiki/Geometric_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param p success probability, \eqn{0 \le} \code{p} \eqn{\le 1}, default value \code{p = 1}.
#' @param type \eqn{standard = 1}, \eqn{shifted = 2}, default \code{type = standard}.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Geometric distribution.
#'
#' @example R/Examples/example_cfN_Geometric.R
#'
#' @export
#'
cfN_Geometric <- function(t, p = 1, type = "standard", cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t)) {
    stop("Enter input parameter t.")
  }

  ## Characteristic function of the (compound) Geometric distribution
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {
    expit <- exp(1i * t)
  } else {
    expit = cfX(t)
  }

  cf <- switch(type,
               standard = 1,
               shifted = expit)

  cf <- cf * (p / (1 - (1 - p) * expit))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
