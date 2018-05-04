#' @title Characteristic function of a linear combination
#' of independent VON MISES random variables
#'
#' @description
#' \code{cf_vonMises(t, mu, kappa, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = sum_{i=1}^N coef_i * X_i} where \eqn{X_i ~ vonMises(\mu_i,\kappa_i)}
#' inedependent RVs, with the locarion parameters \eqn{\mu_i} in Real and the rate
#' parameters \eqn{\kappa_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of the \eqn{vonMises(\mu,\kappa)} distribution is
#' \deqn{cf(t) = cf_vonMises(t,\mu,\kappa) = besseli(t,\kappa)/besseli(0,\kappa) * exp(1i*t*\mu).}
#'
#' @family Continuous Probability Distribution
#' @family Circular Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Von_Mises_distribution}.
#'
#' @param t numerical values (number, vector...).
#' @param mu in \eqn{(-\pi, \pi)}.
#' @param kappa \eqn{> 0}.
#' @param coef vector of the coefficients of the linear combination
#' of the IGamma random variables. If coef is scalar, it is assumed
#' that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y +...+ Y} is
#' sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef_i * X_i}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @details
#' The VON MISES distribution is circular distribution on the interval
#' of length \eqn{2*\pi}, here we consider \eqn{(-\pi,\pi)}, equivalent of the normal
#' distribution with the real parameter \eqn{\mu} and rate parameter \eqn{\kappa > 0}
#' (\eqn{\mu} and \eqn{1/\kappa} are analogous to \eqn{\mu} and \eqn{\sigma^2}, the mean and variance
#' in the normal distribution), on a whole circle, i.e. the interval of angles \eqn{(-\pi,\pi)}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent VON MISES random variables.
#'
#' @example R/Examples/example_cf_vonMises.R
#'
#' @export
#'
cf_vonMises <- function(t,
                        mu = 0,
                        kappa = 1,
                        coef,
                        niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t)) {
    stop("Enter input argument t.")
  }
  if (missing(mu)) {
    mu <- vector()
  }
  if (missing(kappa)) {
    kappa <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- vector()
  }

  ##
  if (length(mu) == 0) {
    mu <- 0
  }
  if (length(kappa) == 0) {
    kappa <- 0
  }
  if (length(coef) == 0) {
    coef <- 1
  }
  if (length(niid) == 0) {
    niid <- 1
  }

  l_max <- max(c(length(mu), length(kappa), length(coef)))
  if (l_max > 1) {
    if (length(mu) == 1) {
      mu <- rep(mu, l_max)
    }
    if (length(kappa) == 1) {
      kappa <- rep(kappa, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(coef, mu, kappa)) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  szt <- dim(t)
  t <- c(t)

  if (length(coef) == 1) {
    aux <- vector()
    for(i in 1:length(t)) {
      aux <- tryCatch(Bessel::BesselI(kappa, abs(t[i] * coef), TRUE) / Bessel::BesselI(kappa, 0, TRUE), error = function(e) 0)
    }
    cf <- aux * exp(1i * t * mu * coef)
  } else {
    aux0 <- matrix(0, dim(t %*% t(coef))[1], dim(t %*% t(coef))[2])
    aux1 <- rep(1, length(t)) %*% t(kappa)
    aux2 <- abs(t %*% t(coef))
    aux3 <- exp(1i * t %*% t(mu * coef))
    for(j in 1:length(t)) {
      for(k in 1:length(coef)) {
        aux0[j,k] <- tryCatch(Bessel::BesselI(aux1[j,k], aux2[j,k], TRUE) / Bessel::BesselI(aux1[j,k], 0, TRUE), error = function(e) 0)
      }
    }
    cf <- apply(aux0 * aux3, 1, prod)
  }

  cf[t == 0] <- 1

  dim(cf) <- szt

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value.")
    }
  }

  return(cf)
}
