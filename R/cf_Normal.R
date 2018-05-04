#' @title Characteristic function of a linear combination of independent NORMAL random variables.
#'
#' @description
#' \code{cf_Normal(t, mu, sigma, coef, niid)} evaluates the characteristic function
#' \eqn{cf(t)} of a linear combination (resp. convolution) of independent NORMAL random variables.
#'
#' That is, \code{cf_Normal} evaluates the characteristic function \eqn{cf(t)}
#' of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ N(\mu_i,\sigma_i)} are inedependent RVs,
#' with means \eqn{\mu_i} and standard deviations \eqn{\sigma_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y} is defined by
#' \deqn{cf(t) = exp(1i*(coef'*\mu)*t - (1/2)*(coef^2'*\sigma^2)*t^2).}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param mu vector of the 'location' parameters \code{mu} in \eqn{R}. If empty, default value is \code{mu = 0}.
#' @param sigma vector of the 'scale' parameters \code{sigma > 0}. If empty, default value is \code{sigma = 1}.
#' @param coef vector of the coefficients of the linear combination of the Normal random variables.
#' If \code{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)} is independently
#' and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @note
#' The characteristic function of a lienar combination of independent
#' Normal random variables has well known exact form.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent NORMAL random variables.
#'
#' @example R/Examples/example_cf_Normal.R
#'
#' @export
#'
cf_Normal <- function(t, mu, sigma, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(mu)) {
    mu <- vector()
  }
  if (missing(sigma)) {
    sigma <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  if (length(sigma) == 0 && length(mu) > 0) {
    sigma <- 1
  } else if (length(sigma) == 0 && length(coef) > 0) {
    sigma <- 1
  } else if (any(sigma == 0) || length(sigma) == 0) {
    sigma <- 1
  }

  if (length(mu) == 0 && length(coef) > 0) {
    mu <- 0
  } else if (length(mu) == 0 && length(sigma) > 0) {
    mu <- 0
  }

  if (length(coef) == 0 && length(sigma) > 0) {
    coef <- 1
  } else if (length(coef) == 0 && length(mu) > 0) {
    coef <- 1
  }

  if (length(niid) == 0) {
    niid <- 1
  }

  l_max <- max(c(length(mu), length(sigma), length(coef)))
  if (l_max > 1) {
    if (length(mu) == 1) {
      mu <- rep(mu, l_max)
    }
    if (length(sigma) == 1) {
      sigma <- rep(sigma, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (any(lengths(list(coef, mu, sigma)) < l_max)) {
      stop("Input size mismatch.")
    }
  }

  szt <- dim(t)
  t <- c(t)

  mean <- as.numeric(t(coef) %*% Conj(mu))
  var <- sum((coef * sigma) ^ 2)
  cf <- exp(1i * mean * t - ((var / 2) * (t ^ 2)))
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
