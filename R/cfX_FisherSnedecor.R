#' @title Characteristic function of the central FISHER-SNEDECOR F-distribution
#'
#' @description
#' \code{cfX_FisherSnedecor(t, df1, df2, coef, niid, tol)} evaluates characteristic function
#' of the central FISHER-SNEDECOR F-distribution with \eqn{df1 > 0} and \eqn{df2 > 0} degrees of freedom.
#'
#' \code{cfX_FisherSnedecor} is an ALIAS NAME of the more general function
#' \code{cf_FisherSnedecor}, used to evaluate the characteristic function
#' of a linear combination of independent FISHER-SNEDECOR F-distributed random variables.
#'
#' The characteristic function of \eqn{X ~ F(df1,df2)} is defined
#' by \eqn{cf(t) = U(df1/2, 1-df2/2, -1i*(df2/df1)*t)},
#' where \eqn{U(a,b,z)} denotes the confluent hypergeometric function of the second kind.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param df1 vector of the  degrees of freedom \code{df1 > 0}. If empty, default value is \code{df1 = 1}.
#' @param df2 vector of the  degrees of freedom \code{df2 > 0}. If empty, default value is \code{df2 = 1}.
#' @param coef vector of the coefficients of the linear combination of the log-transformed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#' @param tol tolerance factor for selecting the Poisson weights, i.e. such that \eqn{PoissProb > tol}.
#' If empty, default value is \code{tol = 1e-12}.
#'
#' @return Characteristic function \eqn{cf(t)} of the central FISHER-SNEDECOR F-distribution.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/F-distribution}.
#'
#' @family Continuous Probability Distribution
#'
#' @references
#' [1] PHILLIPS, P.C.B. The true characteristic function of the F distribution. Biometrika (1982), 261-264.
#'
#' [2] WITKOVSKY, V.: On the exact computation of the density and of the quantiles of linear combinations
#' of t and F random variables. Journal of Statistical Planning and Inference 94 (2001), 1-13.
#'
#' @note Ver.: 06-Oct-2018 17:40:03 (consistent with Matlab CharFunTool v1.3.0, 24-Jun-2017 10:07:43).
#'
#' @example R/Examples/example_cfX_FisherSnedecor.R
#'
#' @export
#'
cfX_FisherSnedecor <- function(t, df1, df2, coef, niid, tol) {

        if(missing(tol)) {
                tol <- numeric()
        }
        if(missing(niid)) {
                niid <- numeric()
        }
        if(missing(coef)) {
                coef <- vector()
        }
        if(missing(df1)) {
                df1 <- vector()
        }
        if(missing(df2)) {
                df2 <- numeric()
        }

        cf <- cf_FisherSnedecor(t, df1, df2, coef, niid, tol)

        return(cf)
}
