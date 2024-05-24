#' @title
#' The confluent hypergeometric function of matrix argument
#'
#' @description
#' \code{Hypergeom1F1Mat(a, b, X, MAX)} computes the confluent hypergeometric function \eqn{1F1(a;b;X)}
#' of a \eqn{(p x p)}-matrix argument \eqn{X}. \code{Hypergeom1F1Mat} is defined for the complex parameters
#' \code{a} and \code{b}, with \eqn{Re(a) > (p-1)/2} and  \eqn{Re(b-a) > (p-1)/2}, and a REAL symmetric matrix argument \code{X}.
#'
#' For more details and definition of the hypergeometric functions with matrix argument see, e.g.,
#' Koev and Edelman (2006) or Muirhead (2009).
#'
#' @references
#' [1] Koev, P. and Edelman, A., 2006. The efficient evaluation of the hypergeometric function of a matrix argument.
#' \emph{Mathematics of Computation}, 75(254), 833-846.
#'
#' [2] Muirhead RJ. Aspects of multivariate statistical theory. John Wiley & Sons; 2009 Sep 25.
#'
#' [3] Butler RW, Wood AT. Laplace approximations for hypergeometric functions with matrix argument.
#' \emph{The Annals of Statistics}. 2002;30(4):1155-77.
#'
#' @param a complex vector of parameters of the hypergeometric function \eqn{1F1^alpha(a;b;X)}.
#' @param b complex vector of parameters  of the hypergeometric function \eqn{1F1^alpha(a;b;X)}.
#' @param X real symmetric \eqn{(p x p)}-matrix argument (alternatively can be specified
#' as a \eqn{(p x p)}-diagonal matrix or a \eqn{p}-vector of the eigenvalues of \eqn{X}).
#' @param MAX maximum number of partitions, \eqn{|\kappa| <= MAX}, default value is \code{MAX = 20}.
#'
#' @family Utility Function
#'
#' @return
#' Hypergeometric sum, \eqn{1F1(a;b;X)}.
#'
#' @note Ver.: 18-Oct-2018 13:16:53 (consistent with Matlab CharFunTool v1.3.0, 25-Oct-2017 14:56:37).
#'
#' @example R/Examples/example_Hypergeom1F1Mat.R
#'
#' @export
#'
Hypergeom1F1Mat <- function(a, b, X, MAX) {

        if(missing(a) || missing(b) || missing(X)) {
                stop("Input parameters a, b, X must be entered.")
        }

        if(missing(MAX)) {
                MAX <- 20
        }

        f <- HypergeompFqMat(a, b, X, alpha = 2, MAX = MAX)

        return(f)
}
