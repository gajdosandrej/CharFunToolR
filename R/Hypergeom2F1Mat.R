#' @title
#' The Gauss hypergeometric function of matrix argument X
#'
#' @description
#' \code{Hypergeom2F1Mat(a, b, c, X, MAX)} computes the The Gauss hypergeometric function \eqn{2F1(a,b;c;X)}
#' of a \eqn{(p x p)}-matrix argument \eqn{X}. \code{Hypergeom2F1Mat} is defined for the complex
#' parameters \code{a}, \code{b}, and \code{c} with \eqn{Re(a) > (p-1)/2} and  \eqn{Re(c-a) > (p-1)/2},
#' and a real symmetric matrix argument \code{X}, with \eqn{Re(X) < I}.
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
#' @param a complex vector of parameters of the hypergeometric function \eqn{2F1(a,b;c;X)}.
#' @param b complex vector of parameters of the hypergeometric function \eqn{2F1(a,b;c;X)}.
#' @param c complex vector of parameters of the hypergeometric function \eqn{2F1(a,b;c;X)}.
#' @param X real symmetric \eqn{(p x p)}-matrix argument (alternatively can be specified
#' as a \eqn{(p x p)}-diagonal matrix or a \eqn{p}-vector of the eigenvalues of \eqn{X}).
#' @param MAX maximum number of partitions, \eqn{|\kappa| <= MAX}, default value is \code{MAX = 20}.
#'
#' @family Utility Function
#'
#' @return
#' Hypergeometric sum, \eqn{2F1(a,b;c;X)}.
#'
#' @note Ver.: 18-Oct-2018 13:52:32 (consistent with Matlab CharFunTool v1.3.0, 22-Aug-2018 12:23:08).
#'
#' @example R/Examples/example_Hypergeom2F1Mat.R
#'
#' @export
#'
Hypergeom2F1Mat <- function(a, b, c, X, MAX) {

        if(missing(a) || missing(b) || missing(c) || missing(X)) {
                stop("Input parameters a, b, c, X must be entered.")
        }

        ## CHECK THE INPUT PARAMETERS
        if(missing(MAX)) {
                MAX <- 20
        }

        sza <- numeric()
        szb <- numeric()
        szc <- numeric()
        if(is.vector(a)) {
                sza <- c(length(a), 1)
        } else {
                sza <- dim(a)
        }
        if(is.vector(b)) {
                szb <- c(length(b), 1)
        } else {
                szb <- dim(b)
        }
        if(is.vector(c)) {
                szc <- c(length(c), 1)
        } else {
                szc <- dim(c)
        }

        l_max <- max(c(length(a), length(b), length(c)))
        if (l_max > 1) {
                if (length(b) == 1) {
                        b <- rep(b, l_max)
                }
                if (length(a) == 1) {
                        a <- rep(a, l_max)
                }
                if (length(c) == 1) {
                        c <- rep(c, l_max)
                }
                if ((any(lengths(list(a, b, c)) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        ## ALGORITHM
        aux_vec <- c(a,b)
        f <- HypergeompFqMat(matrix(aux_vec, length(a), 2), c, X, alpha = 2, MAX = MAX)

        if(!is.null(sza) && max(sza) > 1) {
                dim(f[[1]]) <- sza
                return(f)
        } else if(!is.null(szb) && max(szb) > 1) {
                dim(f[[1]]) <- szb
                return(f)
        } else if(!is.null(szc) && max(szc) > 1) {
                dim(f[[1]]) <- szc
        }

        return(f)
}
