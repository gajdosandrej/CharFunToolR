#' @title Characteristic function of a linear combination
#' of independent VON MISES random variables
#'
#' @description
#' \code{cfC_vonMises(t, mu, kappa, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = sum_{i=1}^N coef_i * X_i} where \eqn{X_i ~ vonMises(\mu_i,\kappa_i)}
#' inedependent RVs, with the locarion parameters \eqn{\mu_i} in Real and the rate
#' parameters \eqn{\kappa_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of the \eqn{vonMises(\mu,\kappa)} distribution is
#' \deqn{cf(t) = cf_vonMises(t,\mu,\kappa) = besseli(t,\kappa)/besseli(0,\kappa) * exp(1i*t*\mu).}
#'
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
#' @example R/Examples/example_cfC_vonMises.R
#'
#' @export
#'
cfC_vonMises <- function(t,
                         mu = 0,
                         kappa = 1,
                         coef,
                         niid) {
  szt <- dim(t)
  t <- c(t)



  cf <-
    unlist(lapply(t, function(t)
      tryCatch((Bessel::BesselI(kappa, abs(t), TRUE) / Bessel::BesselI(kappa, 0, TRUE)) * exp(1i *
                                                                                t * mu),
               error = function(e)
                 0
      )))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
