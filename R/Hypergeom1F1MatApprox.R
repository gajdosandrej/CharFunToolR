#' @title
#'  Computes the approximation of the confluent hypergeometric function \eqn{1F1(a,b,X)} of a matrix argument
#'
#' @description
#' \code{Hypergeom1F1MatApprox(a, b, X)} computes the approximation of the confluent
#' hypergeometric function \eqn{1F1(a,b,X)} of a matrix argument, defined for the complex parameters \code{a} and \code{b},
#' with \eqn{Re(a) > (p-1)/2} and  \eqn{Re(b-a) > (p-1)/2}, and a REAL symmetric \eqn{(p x p)}-matrix argument \code{X}.
#'
#' In fact, \eqn{1F1(a,b,X)} depends only on the eigenvalues of \eqn{X}, so \eqn{X} could be specified
#' as a \eqn{(p x p)}-diagonal matrix or a \eqn{p}-dimensional vector of eigenvalues of the original matrix \eqn{X}, say \eqn{x}.
#'
#' Based on heuristic arguments (not formally proved yet), the value of the confluent
#' hypergeometric function \eqn{1F1(a,b,X)} of a matrix argument is calculated as
#' \eqn{1F1(a;b;X) ~ 1F1(a;b;x(1)) * ... * 1F1(a;b;x(p))},
#' where \eqn{1F1(a;b;x(1))} is the scalar value confluent hypergeometric
#' function \eqn{1F1(a,b,x(i))} with \eqn{[x(1),...,x(p)] = eig(X)}.
#'
#' Here the confluent hypergeometric function \eqn{1F1(a;b;z)} is evaluated for the vector
#' parameters \code{a} and \code{b} and the scalar argument \code{z} by using the simple (4-step) series expansion.
#'
#' @param a complex vector of parameters of the hypergeometric function \eqn{1F1(a;b;X)}.
#' @param b complex vector of parameters  of the hypergeometric function \eqn{1F1(a;b;X)}.
#' @param X real symmetric \eqn{(p x p)}-matrix argument (alternatively can be specified
#' as a \eqn{(p x p)}-diagonal matrix or a \eqn{p}-vector of the eigenvalues of \code{X}), i.e. \eqn{x = eig(X)}.
#'
#' @family Utility Function
#'
#' @return
#' (Approximate) value of the confluent hypergeometric function \eqn{1F1(a;b;X)}, of a matrix argument \code{X}.
#'
#' @note Ver.: 06-Oct-2018 18:45:44 (consistent with Matlab CharFunTool v1.3.0, 19-Jul-2018 16:11:57).
#'
#' @example R/Examples/example_Hypergeom1F1MatApprox.R
#'
#' @export
#'
Hypergeom1F1MatApprox <- function(a, b, X) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(a) || missing(b) || missing(X)) {
                stop("Enter all input parameters.")
        }

        l_max <- max(length(a), length(b))
        if (l_max > 1) {
                if (length(a) == 1) {
                        a <- rep(a, l_max)
                }
                if (length(b) == 1) {
                        b <- rep(b, l_max)
                }
                if ((any(lengths(list(a, b)) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        ## ALGORITHM
        if(is.matrix(X) && dim(X)[1] == dim(X)[2]) {
                x <- eigen(X)
        } else if(is.vector(X)) {
                x <- c(X)
        } else {
                stop("Input size mismatch.")
        }

        p <- length(x)
        f <- 1
        for(i in 1:p) {
                f <- f * Hypergeom1F1SeriesExp(a, b, x[i])
        }

        return(f)
}

Hypergeom1F1SeriesExp <- function(a, b, x, n, tol) {
        # Hypergeom1F1SeriesExp Computes the confluent hypergeometric function
        # 1F1(a;b;z) also known as the Kummer's function M(a,b,z), for the
        # vector parameters a and b and the scalar argument z by using the simple
        # (4-step) series expansion.
        #
        # For more details on confluent hypergeometric function 1F1(a;b;z) or the
        # Kummer's (confluent hypergeometric) function M(a, b, z) see WIKIPEDIA:
        # https://en.wikipedia.org/wiki/Confluent_hypergeometric_function.
        #
        # SYNTAX
        #
        # INPUTS
        #  a      - vector parameter a,
        #  b      - vector parameter b,
        #  x      - scalar argument x,
        #  n      - maximum number of terms used in the series expansion. If empty,
        #           default value of n = 500,
        #  tol    - tolerance for stopping rule. If empty, default value f tol = 1e-14.
        #
        # OUTPUTS
        #  f      - calculated 1F1(a,b,x) of the same dimension as a snd b,
        #  isConv - flag indicator for convergence of the series, isConv = true if
        #           loops < n
        #  loops  - total number of terms used in the series expansion,
        #  n      - used maximum number of terms in the series expansion,
        #  tol    - used tolerance for stopping rule.
        #
        # EXAMPLE 1
        #  t <- seq(-20, 20, length.out = 101)
        #  a <- 1i * t
        #  b <- 2 - 3i * t
        #  x <- 3.14
        #  f <- Hypergeom1F1SeriesExp(a, b, x)

        ## CHECK THE INPUT PARAMETERS
        if(missing(n)) {
                n <- numeric()
        }
        if(missing(tol)) {
                tol <- numeric()
        }

        if(length(n) == 0) {
                n <- 500
        }
        if(length(tol) == 0) {
                tol <- 1e-14
        }

        l_max <- max(c(length(a), length(b)))
        if (l_max > 1) {
                if (length(b) == 1) {
                        b <- rep(b, l_max)
                }
                if (length(a) == 1) {
                        a <- rep(a, l_max)
                }
                if ((any(lengths(list(a, b)) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        f <- 1
        r1 <- 1
        loops <- 0
        for(j in seq(1, n, 4)) {
                loops <- j
                r1 <- r1 * (a + j - 1) / ( j * ( b + j - 1)) * x
                r2 <- r1 * (a + j) / ( (j + 1) * ( b + j)) * x
                r3 <- r2 * (a + j + 1) / ( (j + 2) * ( b + j + 1)) * x
                r4 <- r3 * (a + j + 2) / ( (j + 3) * ( b + j + 2)) * x
                rg <- r1 + r2 + r3 + r4
                f  <- f + rg
                if(max(abs( rg / f )) < tol ) {
                        break
                }

                r1 <-  r4
        }

        isConv <- FALSE
        if(loops < n) {
                isConv <- TRUE
        }

        return(f)
}

