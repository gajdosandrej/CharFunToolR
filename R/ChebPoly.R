#' @title
#' Evaluates the Chebyshev nth order polynomial of the first kind
#'
#' @description
#' \code{ChebPoly(n, x)} evaluates the Chebyshev nth order polynomial of the first kind
#' \eqn{T_n(x)}, where \eqn{T_n(x) = cos(n*acos(x))}.
#'
#' @seealso For more details see also \url{https://en.wikipedia.org/wiki/Chebyshev_polynomials}.
#'
#' @family Utility Function
#'
#' @param n order of polynomial.
#' @param x points, in which we want to compute the value of polynomial.
#'
#' @return  Function returns values of nth order Chebysev polynomial of the first kind evaluated in points \code{x}.
#'
#' @note Ver.: 01-Oct-2018 12:33:49 (consistent with Matlab CharFunTool v1.3.0, 07-Sep-2017 16:20:17).
#'
#' @example R/Examples/example_ChebPoly.R
#'
#' @export
#'
ChebPoly <- function(n, x) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(x)) {
                x <- vector()
        }

        if(length(x) == 0) {
                x <- seq(-1, 1, length.out = 101)
        }

        szx <- dim(x)
        x <- Conj(c(x))

        ## ALGORITHM
        pval  = cos(n %*% t(acos(x)))

        if(length(n) == 1) {
                dim(pval) <- szx
        }
        return(pval)
}
