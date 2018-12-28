#' @title
#' Evaluates the multivariate gamma function of order p
#'
#' @description
#' \code{GammaMulti(z, p)} evaluates the multivariate gamma function of order \code{p},
#' the argument \code{z} may be complex and of any size.
#'
#' @family Utility Function
#'
#' @param z complex argument, of any size (vector, matrix).
#' @param p order of the multivariate gamma, if empty, default value is \code{p = 1}.
#'
#' @return  Function returns values of the multivariate gamma function of order \code{p} evaluated in points \code{z}.
#'
#' @note Ver.: 01-Oct-2018 13:41:36 (consistent with Matlab CharFunTool v1.3.0, 17-Aug-2018 19:45:37).
#'
#' @example R/Examples/example_GammaMulti.R
#'
#' @export
#'
GammaMulti <- function(z, p) {
        ## ALGORITHM

        f <- GammaMultiLog(z, p, 1)

        return(f)
}
