#' @title Characteristic function of a linear combination
#' of independent zero-mean symmetric ARCSINE random variables
#'
#' @description
#' \code{cfS_ArcsineSymmetric(t, coef, niid)} evaluates Characteristic function of a linear combination
#' (resp. convolution) of independent zero-mean symmetric ARCSINE random variables defined on the interval \eqn{(-1,1)}.
#'
#' That is, \code{cfS_ArcsineSymmetric} evaluates the characteristic function
#' \eqn{cf(t) of  Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ ArcsineSymmetric}
#' are independent RVs defined on \eqn{(-1,1)}, for all \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ ArcsineSymmetric} is defined by
#' \deqn{cf(t) = cf_ArcsineSymmetric(t) = besselj(0,t).}
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @importFrom Bessel BesselJ
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic
#' function: An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Arcsine_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param coef vector of the coefficients of the linear combination of the
#' symmetric Arcsinne distributed random variables. If \code{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)} is independently
#' and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent zero-mean symmetric ARCSINE random variables.
#'
#' @example R/Examples/example_cf_ArcsineSymmetric.R
#'
#' @export
#'
cf_ArcsineSymmetric <- function(t, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(coef)) {
    coef <- vector()
  }

  if (missing(niid)) {
    niid <- numeric()
  }

  if (length(coef) == 0) {
    coef <- 1
  }

  if (length(niid) == 0) {
    niid <- 1
  }

  ## Characteristic function
  szt <- dim(t)
  t <- c(t)
  coef <- Conj(t(c(coef)))
  sz_t_coef <- dim(t %*% coef)
  if(length(coef) > 1) {
    t_coef <- t %*% coef
    cf <- matrix(0,sz_t_coef[1],sz_t_coef[2])
    for(i in 1:length(t)) {
      for(j in 1:length(coef)) {
        cf[i,j] <- tryCatch(BesselJ(t_coef[i,j], 0), error = function(e) 0)
      }
    }
    dim(cf) <- sz_t_coef
  } else {
    t_coef <- t %*% coef
    cf <- vector()
    for(k in 1:length(t)) {
      cf <- c(cf, tryCatch(BesselJ(t_coef[k,1], 0), error = function(e) 0))
    }
    dim(cf) <- szt
  }
  # cf <- tryCatch(BesselJ(t %*% coef, 0), error = function(e) 0)
  # cf <- tryCatch(besselJ(t %*% coef, 0), error = function(e) 0)
  # cf <- tryCatch(bessel_Jnu(0, t %*% coef), error = function(e) 0)

  if(length(coef) > 1) {
    cf <- apply(cf, 1, prod)
  }
  cf[t == 0] <- 1

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop('niid should be a scalar (positive integer) value')
    }
  }

  return(cf)

}
