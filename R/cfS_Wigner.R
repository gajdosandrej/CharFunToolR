#' @title Characteristic function of a linear combination (resp. convolution) of
#' the zero-mean symmetric WIGNER distribution defined on the interval \eqn{(-1,1)}.
#'
#' @description
#' \code{cfS_Wigner(t,coef,niid)} is an ALIAS of the more general function \code{cf_WignerSemicircle},
#' used to evaluate the characteristic function of a linear combination of
#' independent symmetric WIGNER SEMICIRCLE distributed random variables.
#'
#' The characteristic function of the symmetric WIGNER distribution on
#' \eqn{(-1,1)} is \eqn{cf(t) = 2*besselj(1,t)/t}.
#'
#' @family Continuous Probability Distribution
#'
#' @references
#'  WITKOVSKY V. (2016). Numerical inversion of a characteristic
#'  function: An alternative tool to form the probability distribution of
#'  output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Wigner_semicircle_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param coef vector of the coefficients of the linear combination of the Beta distributed random variables.
#' If \code{coef} is scalar, it is assumed that all coefficients are equal. If empty, default value is coef = 1.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is sum of \eqn{niid} iid random variables \eqn{Y},
#' where each \eqn{Y = sum_{i=1}^N coef(i) * X_i} is independently and identically distributed random variable.
#' If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function of a linear combination of
#' independent symmetric WIGNER SEMICIRCLE distributed random variables.
#'
#' @note Ver.: 11-Aug-2021 15:51:34 (consistent with Matlab CharFunTool v1.5.1, 18-Sep-2018 00:45:54).
#'
#' @example R/Examples/example_cfS_Wigner.R
#'
#' @export
#
cfS_Wigner <- function(t,coef,niid) {

  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- vector()
  }

  cf<-cf_WignerSemicircle(t,c(),c(),coef,niid)
  return (cf)
}
