#' @title Characteristic function of a symmetric Laplace distribution
#' with scale parameter \eqn{beta>0}.
#'
#' @description
#' \code{cfS_Laplace(t, beta, coef, niid)} evaluates the characteristic function of a linear
#' combination of independent (symmetric) Laplace distributed random variables.
#'
#' That is, \code{cfS_Laplace} evaluates the characteristic function
#' \eqn{cf(t)} of \eqn{Y=sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ Laplace (0,beta_i)}
#' are independent zero-mean RVs with the scale parameters
#' \eqn{beta_i > 0}, for \eqn{i = 1,...,N}
#'
#' @family Continuous Probability Distribution
#' @family Symmetric Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param beta vector of the scale parameters \code{beta > 0}. If empty, default value is \code{beta = 1}.
#' @param coef vector of the coefficients of the linear combination of the LAPLACE random variables.
#' If \code{coef} is scalar, it is assumed that all coefficients are equal. If empty, default value is coef = 1.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is sum of \eqn{niid} iid random variables \eqn{Y},
#' where each \eqn{Y = sum_{i=1}^N coef(i) * X_i} is independently and identically distributed random variable.
#' If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of of a linear
#' combination of independent (symmetric) LAPLACE sistributed random variables.
#'
#' @note Ver.: 10-Aug-2021 17:21:39 (consistent with Matlab CharFunTool v1.5.1, 16-Aug-2018 16:00:43).
#'
#' @example R/Examples/example_cfS_Laplace.R
#'
#' @export
#'

cfS_Laplace <- function(t, beta, coef, niid) {
  if (missing(beta)) {
    beta <- vector()
  }

  if (missing(coef)) {
    coef <- vector()
  }

  if (missing(niid)) {
    niid <- numeric()
  }


  cf<- cf_Laplace(t,c(),beta,coef,niid)


  return(cf)
}
