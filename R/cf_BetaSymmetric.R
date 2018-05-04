#' @title Characteristic function of a linear combination
#' of independent zero-mean symmetric BETA random variables
#'
#' @description
#' \code{cf_BetaSymmetric(t, theta, coef, niid)} evaluates the characteristic function of a linear combination
#' (resp. convolution) of independent zero-mean symmetric BETA random variables defined on the interval  \eqn{(-1,1)}.
#'
#'  That is, \code{cf_BetaSymmetric} evaluates the characteristic function
#'  \eqn{cf(t)} of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i~BetaSymmetric(\theta_i)}
#'  are independent RVs defined on \eqn{(-1,1)}, for all \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ BetaSymmetric(\theta)} is defined by
#' \deqn{cf(t) = cf_BetaSymmetric(t,\theta) = gamma(1/2+\theta) * (t/2)^(1/2-\theta) * besselj(\theta-1/2,t).}
#'
#' SPECIAL CASES: \cr
#' 1) \eqn{\theta = 1/2}; Arcsine distribution on \eqn{(-1,1)}:     \eqn{cf(t) = besselj(0,t)}, \cr
#' 2) \eqn{\theta = 1};   Rectangular distribution on \eqn{(-1,1)}: \eqn{cf(t) = sin(t)/t}.
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
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param theta the 'shape' parameter \code{theta > 0}. If empty, default value is \code{theta = 1}.
#' @param coef vector of the coefficients of the linear combination of Beta distributed random variables. If \code{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)} is independently
#' and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent zero-mean symmetric BETA random variables.
#'
#' @example R/Examples/example_cf_BetaSymmetric.R
#'
#' @export
#'
cf_BetaSymmetric <- function(t, theta, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(theta)) {
    theta <- vector()
  }

  if (missing(coef)) {
    coef <- vector()
  }

  if (missing(niid)) {
    niid <- numeric()
  }

  if (length(theta) == 0) {
    theta <- 1
  }

  if (length(coef) == 0) {
    coef <- 1
  }

  if (length(niid) == 0) {
    niid <- 1
  }

  l_max <- max(c(length(theta), length(coef)))
  if (l_max > 1) {
    if (length(theta) == 1) {
      theta <- rep(theta, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(coef, theta)) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function
  szt <- dim(t)
  t <- c(t)
  cf <- vector()

  if (length(coef) == 1) {
    # cf <- vector()
    coef_t <- coef * t
    for(i in 1:length(t)) {
      cf <- c(cf, tryCatch(Re(exp(GammaLog(0.5 + theta) + (0.5 - theta) * log(0.5 * coef_t[i] + 0i)) * Bessel::BesselJ(coef_t[i], theta - 0.5)), error = function(e) 0))
    }
  } else {
    aux1 <- t %*% t(coef)
    aux2 <- rep(1, length(t)) %*% t(theta)
    cf <- matrix(0, dim(aux1)[1], dim(aux1)[2])
    for(j in 1:length(t)) {
      for(k in 1:length(coef)) {
        cf[j,k] <- tryCatch(Re(exp(GammaLog(0.5 + aux2[j,k]) + (0.5 - aux2[j,k]) * log(0.5 * aux1[j,k] + 0i)) * Bessel::BesselJ(aux1[j,k], aux2[j,k] -0.5)), error = function(e) 0)
      }
    }
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
