#' @title Characteristic function of a linear combination
#' of independent zero-mean symmetric RECTANGULAR random variables
#'
#' @description
#' \code{cf_RectangularSymmetric(t, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of a linear combination (resp. convolution)of independent zero-mean symmetric RECTANGULAR
#' random variables defined on the interval \eqn{(-1,1)}.
#'
#' That is, \code{cf_RectangularSymmetric} evaluates the characteristic function
#' \eqn{cf(t)} of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ RectangularSymmetric}
#' are independent uniformly distributed RVs defined on \eqn{(-1,1)}, for all \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ RectangularSymmetric} is defined by
#' \deqn{cf(t) = cf_RectangularSymmetric(t) = sinc(t) = sin(t)/t.}
#'
#' @family Continuous Probability Distribution
#' @family Symmetric Probability Distribution
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic
#' function: An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param coef vector of the coefficients of the linear combination
#' of the Beta distributed random variables. If coef is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is
#' sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#'  is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent zero-mean symmetric RECTANGULAR random variables.
#'
#' @example R/Examples/example_cf_RectangularSymmetric.R
#'
#' @export
#'
cf_RectangularSymmetric <- function(t, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }
  if (length(coef) == 0) {
    coef <- 1
  }
  if (length(niid) == 0) {
    niid <- 1
  }

  ## Characteristic function
  szt <- dim(t)
  t <- c(t)
  coef <- t(Conj(coef))
  cf <- apply(sin(t %*% coef) / (t %*% coef), 1, prod)

  dim(cf) <- szt
  cf[t == 0] <- 1

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value.")
    }
  }

  return(cf)

}
