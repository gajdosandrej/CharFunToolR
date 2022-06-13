#' @title
#' Evaluates function by using Chebyshev polynomial approximation
#'
#' @description
#'   \code{InterpChebValues(values,x,domain)} evaluates the function fun by using Chebyshev
#'   polynomial approximation specified by the function values evaluated at
#'   Chebyshev points with specified domain.
#'
#' @family Utility Function
#'
#' @param values function values evaluated at Chebyshev points of the 2nd kind.
#' @param x points, in which we want to compute the value of polynomial.
#' @param domain vector containing lower and upper bound of interval.
#'
#' @return Function returns the function fun by using Chebyshev
#'          polynomial approximation specified by the function values evaluated at
#'          Chebyshev points with specified domain.
#'
#' @note Ver.: 16-Nov-2021 15:52:07 (consistent with Matlab CharFunTool v1.3.0, 28-May-2021 14:28:24).
#'
#' @example R/Examples/example_ChebPolyValues.R
#'
#' @export

InterpChebValues <- function(values,x,domain){
        if(missing(domain)) {
                domain <- vector()
        }
        if(missing(x)) {
                x <- vector()
        }
        coeffs <- ChebCoefficients(values)
        fun_x_domain <- ChebPolyValues(coeffs,x,domain)

        return(fun_x_domain)

}
