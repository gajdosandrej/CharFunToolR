#' @title Characteristic function of a linear combination (resp. convolution) of
#' independent WIGNER SEMICIRCLE random variables  defined on the interval
#' \eqn{(mu-R,mu+R)}.
#'
#' @description
#'  That is, \code{cf_WignerSemicircle(t,mu,R,coef,niid)} evaluates the characteristic function
#'  \eqn{cf(t)} of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ WignerSemicircle}
#'  are independent RVs defined on \eqn{(mu_i-R_i,mu_i+R_i)}, for all \eqn{i = 1,...,N}.
#'
#'  The characteristic function of \eqn{X ~ WignerSemicircle(mu,R)} is defined by
#'  \eqn{cf(t) = 2*exp(1i*t*mu).*besselj(1,R*t)/(R*t)};
#'
#'  @family Continuous Probability Distribution
#'
#'  @seealso For more details see WIKIPEDIA:
#'  \url{https://en.wikipedia.org/wiki/Wigner_semicircle_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param mu vector of the 'location' parameters mu in R. If empty, default value is \code{mu = 0}.
#' @param R vector of the 'radius' prameters \eqn{R > 0}. If empty, default value is \code{R =1}.
#' @param coef vector of the coefficients of the linear combination of the Beta distributed random variables.
#' If \code{coef} is scalar, it is assumed that all coefficients are equal. If empty, default value is coef = 1.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is sum of \eqn{niid} iid random variables \eqn{Y},
#' where each \eqn{Y = sum_{i=1}^N coef(i) * X_i} is independently and identically distributed random variable.
#' If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear
#' combination (res. convolution) of independent WIGNER SEMICIRCLE random variables.
#'
#' @note Ver.: 10-Aug-2021 17:44:17 (consistent with Matlab CharFunTool v1.5.1, 18-Sep-2018 00:32:58).
#'
#' @example R/Examples/example_cf_WignerSemicircle.R
#'
#' @export
#'
 cf_WignerSemicircle <- function(t,mu,R,coef,niid) {

  ## CHECK THE INPUT PARAMETERS
  if(missing(mu)) {
    mu <- vector()
  }
  if(missing(R)) {
    R <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- numeric()
  }

  ##
  if(length(mu) == 0 ) {
    mu <- 0
  }
  if(length(R) == 0 ) {
    R <- 1
  }
  if(length(coef) == 0 ){
    coef <- 1
  }
  if(length(niid) == 0 ){
    niid <- 1
  }

  ## Characteristic function

  l_max <- max(c(length(coef), length(mu), length(R)))
  if (l_max > 1) {
    if (length(mu) == 1) {
      mu <- rep(mu, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (length(R) == 1) {
      R <- rep(R, l_max)
    }
    if ((any(lengths(list(coef, mu, R)) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  szt <- dim(t)
  t <- c(t)
  t_mu_coef <- t%*%t(mu*coef)
  t_R_coef <- t%*%t(R*coef)

  cf<-tryCatch(2*exp(1i*t_mu_coef)* (Bessel::BesselJ(t_R_coef, 1))/(t_R_coef), error = function(e) 0)

  cf<-apply(cf, 1, prod)

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





