#' @title Characteristic function of a linear combination
#' of independent LOG-TRANSFORMED WILK's LAMBDA distributed random variables
#'
#' @description
#' \code{cf_LogRV_WilksLambdaNC(t, p, n, q, delta, coef, MAX)} evaluates
#' characteristic function of a linear combination (resp. convolution)
#' of independent LOG-TRANSFORMED WILK's LAMBDA distributed random variables,
#' with non-central distributions specified by the parameters \code{p, n, q},
#' the non-centrality parameters \code{delta}, and the coefficients \code{coef}.
#'
#' That is, \code{cf_LogRV_WilksLambdaNC} evaluates the characteristic function
#' of a random variable \eqn{Y  = coef_1*W_1 +...+ coef_N*W_N, such that cf_Y(t) = cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t)},
#' where \eqn{cf_W_i(t)} is CF of \eqn{W_i = log(Lambda_i)}, and each Lambda_i has non-central WILK's LAMBDA distribution,
#' \eqn{Lambda_i ~ Lambda(p_i,m_i,n_i,delta_i)}, for \eqn{i = 1,...,N}.
#'
#' In particular, \eqn{\Lambda_i = det(E_i)/det(E_i + H_i)}, with  \eqn{E_i} and \eqn{H _i}
#' being independent random matrices with Wishart distributions, \eqn{E_i ~ Wp_i(m_i,\Sigma_i)}
#' with central Wishart distribution, and \eqn{H_i ~ Wp_i(n_i,\Sigma_i,delta)} with non-central Wishart distribution,
#' with \eqn{n_i >= p_i} and  \eqn{\Sigma_i > 0} (an unknown positive definite symmetric covariance matrix), for all \eqn{i = 1,...,N}.
#'
#' The non-central distribution of MINUS LOG-TRANSFORMED WILK's LAMBDA STATISTIC,
#' say \eqn{L ~ Lambda(p,n,q,delta)}, with \eqn{L} in \eqn{(0,1)}, is specified by its characteristic function
#' \deqn{cf_{-log(L)}(t) = cf_LogRV_Beta(-t,n/2,q/2) * ... * cf_LogRV_Beta(-t,(n+1-p)/2,q/2) * Hypergeom1F1Mat(-i*t,-i*t+(n+q)/2,-delta/2)},
#' i.e. the characteristic function of the non-central distribution differs
#' from the central distribution by a factor specified by the confluent
#' hypergeometric function \eqn{1F1(a;b;X)} of a matrix argument \eqn{X},
#' with the complex (vector) parameters specified as \eqn{a = -i*t, b = -i*t+(n+q)/2},
#' and the matrix argument \eqn{X = -delta/2}, where \eqn{delta} is the non-centrality (matrix) parameter of the distribution.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param p vector of the dimension parameters \eqn{p = (p_1,...,p_N)}.
#' If empty, default value is \code{p = (1,...,1)}.
#' @param n vector of degrees of freedom of the Wishart matrices \eqn{E_i},
#' \eqn{n = (n_1,...,n_N)}. If empty, default value is \code{n = (1,...,1)}.
#' @param q vector of degrees of freedom of the Wishart matrices \eqn{H_i},
#' \eqn{q = (q_1,...,q_N)}. If empty, default value is \code{q = (1,...,1)}.
#' @param delta p.s.d. matrix or vector of nonnegative eigenvalues or
#' array of matrices or vectors of eigenvalues of the non-centrality parameters,
#' \eqn{delta = (delta_1,...,delta_N)}. Default value is an empty array.
#' @param coef vector of the coefficients of the linear combination of the log-transformed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal.
#' If empty, default value is \code{coef = 1}.
#' @param MAX the maximum number of partitions used for computing
#' the hypergeometric \eqn{1F1} function with matrix argument, for more details see \code{HypergeompFqMat}.
#' If MAX is empty, default value is set to \code{MAX = 0}. If \code{MAX = 0},
#' then the algorithm uses the fast (approximate) method for computing the confluent
#' hypergeometric function \eqn{1F1(a;b;X)}, see \code{Hypergeom1F1MatApprox},
#' otherwise the algorithm uses the algorithm \code{Hypergeom1F1Mat}
#' based on computing the truncated expansion series of \eqn{1F1(a;b;X)} with \eqn{MAX} number of partitions.
#'
#' @details
#' Computing CF of the LOG-TRANSFORMED NON-CENTARL WILK's LAMBDA random
#' variable depends on computing the generalized hypergeometric function of matrix argument.
#' \code{cf_LogRV_WilksLambdaNC} uses modified implementation for computing
#' the truncated hypergeometric function, see the algorithm \code{HypergeompFqMat},
#' as originaly suggested in Koev and Edelman (2006). The truncation
#' of the hypergeometric series is controled by the parameter MAX (start with \eqn{MAX = 20}).
#' If necessary, the parameter MAX should be set to sufficiently large value but this can be numerically intractable.
#' If \eqn{MAX = 0} (now it is the default value), the algorithm uses the fast
#' (approximate) method for computing the confluent hypergeometric function
#' \eqn{1F1(a;b;X)}, see \code{Hypergeom1F1MatApprox}, otherwise the algorithm
#' uses the algorithm \code{Hypergeom1F1Mat} based on computing the truncated expansion
#' series of \eqn{1F1(a;b;X)} with \eqn{MAX} number of partitions.
#' The approximate alternative method (if \eqn{MAX = 0}) is based on heuristic
#' arguments (not formally proven), and the value of the confluent
#' hypergeometric function \eqn{1F1(a,b,X)} of a matrix argument
#' is calculated (approximately) as \eqn{1F1(a;b;X) ~ 1F1(a;b;x(1)) * ... * 1F1(a;b;x(p))},
#' where \eqn{1F1(a;b;x(1))} is the scalar value confluent hypergeometric
#' function \eqn{1F1(a,b,x(i))} with \eqn{(x[1],...,x[p]) = eigen(X)}.
#'
#' By default (or if \code{MAX} is set to value \code{MAX = 0}), \code{cf_LogRV_WilksLambdaNC}
#' uses this alternative method for computing \eqn{1F1(a;b;X)}.
#'
#' @importFrom plyr alply
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent
#' LOG-TRANSFORMED WILK's LAMBDA distributed random variables.
#'
#' @references
#' [1] Koev, P. and Edelman, A., 2006. The efficient evaluation
#' of the hypergeometric function of a matrix argument.
#' Mathematics of Computation, 75(254), 833-846.
#'
#' [2] Witkovsky, V., 2018. Exact distribution of selected multivariate test
#' criteria by numerical inversion of their characteristic functions.
#' arXiv preprint arXiv:1801.02248.
#'
#' @seealso For more details see WIKIPEDIA:
#' https://en.wikipedia.org/wiki/Wilks%27s_lambda_distribution.
#'
#' @family Continuous Probability Distribution
#'
#' @note Ver.: 05-Oct-2018 11:46:29 (consistent with Matlab CharFunTool v1.3.0, 10-Aug-2018 15:46:49).
#'
#' @example R/Examples/example_cf_LogRV_WilksLambdaNC.R
#'
#' @export
#'
cf_LogRV_WilksLambdaNC <- function(t, p, n, q, delta, coef, MAX) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(MAX)) {
                MAX <- numeric()
        }
        if(missing(coef)) {
                coef <- vector()
        }
        if(missing(delta)) {
                delta <- list()
        }
        if(missing(q)) {
                q <- vector()
        }
        if(missing(n)) {
                n <- vector()
        }
        if(missing(p)) {
                p <- vector()
        }

        if(length(p) == 0) {
                p <- 1
        }
        if(length(n) == 0) {
                n <- 1
        }
        if(length(q) == 0) {
                q <- 1
        }
        # if(length(delta) == 0) {
        #         delta <- list()
        # }
        if(length(coef) == 0) {
                coef <- 1
        }
        if(length(MAX) == 0) {
                MAX <- 0
        }

        if(!is.list(delta)) {
                dim_delta <- dim(delta)
                if(!is.null(dim_delta) && dim_delta[1] == dim_delta[2]) {
                        delta <- list(eigen((delta + t(Conj(delta))) / 2))
                } else if (is.vector(delta)) {
                        delta <- list(delta)
                } else {
                        stop('Input size mismatch.')
                }
        } else {
                if(length(delta) > 0) {
                        for(i in 1:length(delta)) {
                                dim_delta <- dim(delta[[i]])
                                if(!is.null(dim_delta) && dim_delta[1] == dim_delta[2] && dim_delta[1] > 0) {
                                        delta[[i]] <- eigen((delta[[i]] + t(Conj(delta[[i]]))) / 2)
                                } else if (!is.null(dim_delta) && min(dim_delta[1], dim_delta[2]) > 1) {
                                        stop('Input size Mismatch.')
                                } # one more else if for delta containing vector ?
                        }
                }
        }

        # if(is.matrix(delta) && is.list(delta)) {
        #         for(i in 1:ncol(delta)) {
        #                 dim_deltai <- dim(delta[[1,i]])
        #                 if(!is.null(dim_deltai) && dim_deltai[1] == dim_deltai[2] && dim_deltai[1] > 0) {
        #                         delta[[1,i]] <- eigen((delta[[1,i]] + t(Conj(delta[[1,i]]))) / 2)
        #                 } else if(!is.null(dim_deltai) && min(dim_deltai[1], dim_deltai[2]) > 1) {
        #                         stop('Input size Mismatch.')
        #                 }
        #         }
        # } else {
        #         if(is.matrix(delta)) {
        #                 dim_delta <- dim(delta)
        #         }
        # }

        ## Check size of the parameters
        l_max <- max(c(length(p), length(n), length(q), length(coef)))
        if (l_max > 1) {
                if (length(p) == 1) {
                        p <- rep(p, l_max)
                }
                if (length(n) == 1) {
                        n <- rep(n, l_max)
                }
                if (length(q) == 1) {
                        q <- rep(q, l_max)
                }
                if (length(coef) == 1) {
                        coef <- rep(coef, l_max)
                }
                if ((any(lengths(list(
                        coef, q, n, p
                )) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        # if(length(delta) == 1) {
        #         delta <- plyr::alply(replicate(length(coef), delta[[1]]), 3)
        # }

        if(length(delta) == 1) {
                if(length(coef) > 1) {
                        for(i in 2:length(coef)) {
                                delta[[i]] <- delta[[1]]
                        }
                }
        }

        ## Characteristic function of a linear combination
        szt <- dim(t)
        t <- c(t)

        cf <- 1
        for(i in 1:length(coef)) {
                alpha <- (n[i] + 1 - (1:p[i])) / 2
                beta <- q[i] / 2
                cf <- cf * cf_LogRV_Beta(coef[i] * t, alpha, beta)
                if(length(delta) > 0) {
                        if(length(delta[[i]]) > 0) {
                                if(MAX == 0) {
                                        if(max(t) > 1e298) {
                                                cf <- NaN
                                        } else {
                                                cf <- cf * Hypergeom1F1MatApprox(1i*coef[i]*t,
                                                                                 1i*coef[i]*t + (n[i]+q[i])/2, -delta[[i]]/2)
                                        }
                                } else if (MAX > 0) {
                                        cf <- cf * Hypergeom1F1Mat(1i*coef[i]*t,
                                                                   1i*coef[i]*t + (n[i]+q[i])/2, -delta[[i]]/2,ceiling(MAX))[[1]]
                                } else {
                                        cf <- cf * Hypergeom1F1Mat(1i*coef[i]*t,
                                                                   1i*coef[i]*t + (n[i]+q[i])/2, -delta[[i]]/2,0)[[1]]
                                }
                        }
                }
        }

        dim(cf) <- szt
        cf[t==0] <- 1

        return(cf)
}
