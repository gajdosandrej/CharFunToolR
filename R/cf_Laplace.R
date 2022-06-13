#' @title  Characteristic function of a linear combination (resp. convolution) of
#' independent LAPLACE random variables with location parameter \eqn{mu} (real)
#' and scale parameter \code{beta > 0}.
#'
#' @description \code{cf_Laplace(t,mu,beta,coef,niid)} evaluates the characteristic function
#' \eqn{cf(t)} of \eqn{Y = \sum_{i=1}^N coef_i * X_i}, where  \eqn{X_i ~ Laplace (\mu_i,\beta_{i})}
#' are inedependent RVs, with real location parameters \eqn{\mu_{i}} and the scale parameters \eqn{\beta_{i} > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y} is defined by
#' \deqn{ cf(t)=Prod(exp(li* t * coef(i) \eqn {mu(i)} ) / (1+(t*coef(i)*\eqn{beta(i)})^2) )}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Laplace_distribution}.
#'

#'@param t vector or array of real values, where the CF is evaluated.
#'@param mu vector of real location parameters. If empty, default value is \code{mu = 0}.
#'@param beta vector of the scale parameters \code{beta > 0}. If empty, default value is \code{beta = 1}.
#'@param coef vector of the coefficients of the linear combination of the LAPLACE random variables.
#'If \code{coef} is scalar, it is assumed that all coefficients are equal. If empty, default value is coef = 1.
#'@param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is sum of \eqn{niid} iid random variables \eqn{Y},
#'where each \eqn{Y = sum_{i=1}^N coef(i) * X_i} is independently and identically distributed random variable.
#'If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent LAPLACE random variables.
#'
#' @note Ver.: 08-Aug-2021 16:19:30 (consistent with Matlab CharFunTool v1.5.1, 16-Aug-2018 16:00:43).
#'
#' @example R/Examples/example_cf_Laplace.R
#'
#' @export
#'
cf_Laplace <- function(t,mu,beta,coef,niid) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(mu)) {
    mu <- vector()
  }
  if(missing(beta)) {
    beta <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- vector()
  }

  ##
    if(length(mu) == 0 ) {
      mu <- 0
    }
    if(length(beta) == 0 ) {
      beta <- 1
    }
    if(length(coef) == 0 ){
      coef <- 1
    }


  ## Equal size of the parameters

  l_max <- max(c(length(coef), length(mu), length(beta)))
  if (l_max > 1) {
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (length(mu) == 1) {
      mu <- rep(mu, l_max)
    }
    if (length(beta) == 1) {
      beta <- rep(beta, l_max)
    }
    if ((any(lengths(list(coef, mu, beta)) < l_max))) {
      stop("Input size mismatch.")
    }
  }


## Characteristic function
  szt <- dim(t)
  t <- c(t)

   coef_mu <- coef*mu
   coef_beta <- coef*beta
#cf<-exp(t(apply(as.matrix(1i*t), 2, FUN = '*', coef_mu)))/(1+t(apply(as.matrix(t^2), 2, FUN = '*', coef_beta^2)))
cf<-(exp(1i*t%*%t(coef_mu)))/(1+(t^2%*%t(coef_beta^2)))

cf <- apply(cf, 1, prod)

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

