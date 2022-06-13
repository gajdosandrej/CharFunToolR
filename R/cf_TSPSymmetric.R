#' @title  Characteristic function of a linear combination (resp. convolution) of
#' independent symmetric (location and scale shifted) TWO-SIDED-POWER (TSP)
#' random variables.
#'
#' @description
#' That is, \code{cfS_TSPSymmetric(t, theta, mu, sigma, coef, niid)} evaluates the characteristic function \eqn{cf(t)} of
#' \eqn{Y = sum_{i=1}^N coef_i * (mu_i + sigma_i * X_i)}, where \eqn{X_i ~
#' TSP(theta_i)} are inedependent RVs, with symmetric TSP distributions
#' defined on the interval \eqn{(-1,1)} with zero mean and variance \eqn{Var(X_i) =
#' 2*theta_i*gamma(theta_i)/gamma(3+theta_i)}, where \code{theta_i > 0} are shape
#' parameters for all zeqn{i = 1,...,N}.
#'
#'  The characteristic function of the random variable\eqn{mu + sigma*X}, where
#'   X ~ TSP(theta) is given by
#'   \code{cf(t) = cfS_TSPSymmetric(t,theta,mu,sigma) = 1/2 * exp(1i*t*mu) * (hypergeom1F1(1,1+theta,1i*t*sigma) + hypergeom1F1(1,1+theta,-1i*t*sigma))}.
#'
#'   Hence, the characteristic function of \eqn{Y  = coef_1*(mu_1+sigma_1*X_1) + coef_N*(mu_N+sigma_N*X_N)} is \eqn{cf_Y(t) = exp(1i*mu*t) *(cf_1(coef_1*sigma_1*t) * cf_N(coef_N*sigma_N*t))}, where \eqn{cf_i(t)} is
#'   the characteristic function of \eqn{X_i ~ TSP(theta_i)}.
#'
#' SPECIAL CASES: \cr
#' 1) \eqn{\theta = 1/2}; Arcsine distribution on \eqn{(-1,1)}:     \eqn{cf(t) = besselj(0,t)}, \cr
#' 2) \eqn{\theta = 1};   Rectangular distribution on \eqn{(-1,1)}: \eqn{cf(t) = sin(t)/t}.
#'
#' @family Continuous Probability Distribution
#' @family Symmetric Probability Distribution
#'
#' @references
#'  [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
#'  function: An alternative tool to form the probability distribution of
#'  output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'  [2] VAN DORP, R.J., KOTZ, S. (2003). Generalizations of two-sided power
#'  distributions and their convolution. Communications in
#'  Statistics-Theory and Methods, 32(9), 1703-1723.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param theta vector of the shape parameters \code{theta > 0}. If theta is scalar, it is assumed that all parameters theta are equal.
#' If empty, default value is \code{theta = 1}.
#' @param mu  vector of location parameters, mu in Real. If empty, default value is \code{mu = 0}.
#' @param sigma vector of scale parameters, \code{sigma_i > 0}. If empty, default value is \code{sigma = 1}.
#' @param coef  vector of the coefficients of the linear combination of the log-transformed random variables. If coef is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient n, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)} is independently
#' and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function of a linear combination (resp. convolution) of independent symmetric (location and scale shifted) TWO-SIDED-POWER (TSP) random variables.
#'
#' @note Ver.: 11-Aug-2021 16:18:40 (consistent with Matlab CharFunTool v1.5.1, 24-Jun-2017 18:25:56).
#'
#' @example R/Examples/example_cf_TSPSymmetric.R
#'
#' @export
#





cf_TSPSymmetric <- function(t, theta, mu, sigma, coef, niid) {

## CHECK THE INPUT PARAMETERS
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
  ##
    if(length(theta) == 0 ) {
    theta <- 1
  }
    if(length(mu) == 0) {
    mu <- 0
  }
  if(length(sigma) == 0){
    sigma <- 1
  }
    if(length(coef) == 0){
    coef <- 1
  }
  if(length(niid) == 0){
    niid <- 1
  }

  ## Equal size of the parameters

  l_max <- max(c(length(coef), length(theta), length(mu), length(sigma)))
  if (l_max > 1) {
    if (length(mu) == 1) {
      mu <- rep(mu, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (length(theta) == 1) {
      theta <- rep(theta, l_max)
    }
    if (length(sigma) == 1) {
      sigma <- rep(sigma, l_max)
    }
    if ((any(lengths(list(coef, mu, theta, sigma)) < l_max))) {
      stop("Input size mismatch.")
    }
  }
  ## Characteristic function
  szt <- dim(t)
  t <- abs(c(t))

  cf<-1


  for(i in 1:length(coef)) {

  cf <-cf *exp(1i*t*mu[i])*(hypergeom1F1(1, 1+theta[i],1i*t*coef[i]*sigma[i])$f+hypergeom1F1(1,1+theta[i],-1i*t*coef[i]*sigma[i])$f)/2
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


