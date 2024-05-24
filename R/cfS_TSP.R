#' @title Characteristic function of the TWO-SIDED-POWER (TSP) distribution with shape parameter \code{theta > 0}.
#'
#' @description
#' \code{cfS_TSP(t, theta, mu, sigma, coef, niid)} is an ALIAS of the more general function \code{cf_TSPSymmetric}, used to evaluate the characteristic function of a linear combination of
#' independent (location-scale) TSP distributed random variables.
#'
#'  The characteristic function of the random variable \eqn{X ~ TSP(theta)} is \eqn{cf(t) = 1/2 * (hypergeom1F1(1,1+theta,1i*t) + hypergeom1F1(1,1+theta,-1i*t))}.
#'
#' #' SPECIAL CASES: \cr
#' 1) \eqn{\theta = 1/2}; Arcsine distribution on \eqn{(-1,1)}:     \eqn{cf(t) = besselj(0,t)}, \cr
#' 2) \eqn{\theta = 1};   Rectangular distribution on \eqn{(-1,1)}: \eqn{cf(t) = sin(t)/t}.
#'
#' @family Continuous Probability Distribution
#'
#' @references
#' [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
#' function: An alternative tool to form the probability distribution of
#' output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#' [2] VAN DORP, R.J., KOTZ, S. (2003). Generalizations of two-sided power
#' distributions and their convolution. Communications in
#' Statistics-Theory and Methods, 32(9), 1703-1723.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param theta the shape parameter, \code{theta > 0}. If empty, default value is \code{theta = 1}.
#' @param mu  vector of location parameters, mu in Real. If empty, default value is \code{mu = 0}.
#' @param sigma vector of scale parameters, \code{sigma_i > 0}. If empty, default value is \code{sigma = 1}.
#' @param coef  vector of the coefficients of the linear combination of the log-transformed random variables. If \eqn{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient n, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)} is independently
#' and identically distributed random variable. If empty, default value is \code{niid = 1}.

#' @return Characteristic function of a linear combination  of independent (location scale TSP distributed random variables.
#'
#' @note Ver.: 11-Aug-2021 16:37:22 (consistent with Matlab CharFunTool v1.5.1, 24-Jun-2017 18:25:56).
#'
#' @example R/Examples/example_cfS_TSP.R
#'
#' @export
#
cfS_TSP <- function(t, theta, mu, sigma, coef, niid){

  if(missing(theta)) {
    theta <- vector()
  }
  if(missing(mu)) {
    mu <- vector()
  }
  if(missing(sigma)) {
    sigma <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- numeric()
  }

  cf<-cf_TSPSymmetric(t,theta,mu,sigma,coef,niid)
  return(cf)
}
