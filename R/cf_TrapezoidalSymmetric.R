#' @title Characteristic function of a linear combination
#' of independent zero-mean symmetric TRAPEZOIDAL random variables
#'
#' @description
#' \code{cf_TrapezoidalSymmetric(t, lambda, coef, niid)} evaluates
#' the characteristic function of a linear combination (resp. convolution)
#' of independent zero-mean symmetric TRAPEZOIDAL random variables defined on the interval \eqn{(-1,1)}.
#'
#' That is, \code{cf_TrapezoidalSymmetric} evaluates the characteristic function
#' \eqn{cf(t)} of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ TrapezoidalSymmetric(\lambda_i)}
#' are independent RVs defined on \eqn{(-1,1)}, for all \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ TrapezoidalSymmetric(\lambda)}, with
#' \eqn{E(X) = 0} and \eqn{Var(X) = (1+\lambda^2)/6}, is
#' \deqn{cf(t) = cf_TrapezoidalSymmetric(t) = cf_RectangularSymmetric(w*t))*cf_RectangularSymmetric((1-w)*t) = (sin(w*t)/(w*t))*(sin((1-w)*t)/((1-w)*t)),}
#' where \eqn{w  = (1+\lambda)/2}.
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Trapezoidal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param lambda parameter of the offset, \eqn{0 \le} \code{lambda} \eqn{\le 1}.
#' If empty, default value is \code{lambda = 0}
#' @param coef vector of the coefficients of the linear combination
#' of the zero-mean symmetric TRAPEZOIDAL random variables. If coef is scalar,
#' it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent zero-mean symmetric TRAPEZOIDAL random variables.
#'
#' @example R/Examples/example_cf_TrapezoidalSymmetric.R
#'
#' @export
#'
cf_TrapezoidalSymmetric <- function(t, lambda, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(lambda)) {
    lambda <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  if (length(lambda) == 0) {
    lambda <- 1
  }
  if (length(coef) == 0) {
    coef <- 1
  }
  if (length(niid) == 0) {
    niid <- 1
  }

  l_max <- max(c(length(lambda), length(coef)))
  if (l_max > 1) {
    if (length(lambda) == 1) {
      lambda <- rep(lambda, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(coef, lambda)) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function
  szt <- dim(t)
  t <- c(t)
  w1 <- coef * (1 / 2 + lambda / 2)
  w2 <- coef * (1 / 2 - lambda / 2)

  cf <- apply((sin(t %*% t(w1)) / (t %*% t(w1))) * (sin(t %*% t(w2)) / (t %*% t(w2))), 1, prod)
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
