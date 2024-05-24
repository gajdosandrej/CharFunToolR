#' @title Characteristic function of of a linear combination
#' of independent zero-mean symmetric TRIANGULAR random variables
#'
#' @description
#' \code{cf_TriangularSymmetric(t, coef, niid)} evaluates the characteristic function
#' of a linear combination (resp. convolution) of independent zero-mean symmetric TRIANGULAR
#' random variables defined on the interval \eqn{(-1,1)}.
#'
#' That is, \code{cf_TriangularSymmetric} evaluates the characteristic function
#' \eqn{cf(t)} of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ TriangularSymmetric}
#' are independent uniformly distributed RVs defined on \eqn{(-1,1)}, for all \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ TriangularSymmetric} is defined
#' by \deqn{cf(t) = cf_TriangularSymmetric(t) = (2-2*cos(t))/t^2.}
#'
#' @family Continuous Probability Distribution
#' @family Symmetric Probability Distribution
#'
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Triangular_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param coef vector of the coefficients of the linear combination
#' of the zero-mean symmetric TRIANGULAR random variables. If coef is scalar,
#' it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent zero-mean symmetric TRIANGULAR random variables.
#'
#' @note Ver.: 16-Sep-2018 18:40:04 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cf_TriangularSymmetric.R
#'
#' @export
#'
cf_TriangularSymmetric <- function(t, coef, niid) {
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
  coef <- Conj(coef)
  cf <- (2 - 2 * cos(t %*% t(coef))) / ((t %*% t(coef)) ^ 2)
  if(length(coef) > 1){
    cf <- apply(cf, 1, prod)
  }
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
