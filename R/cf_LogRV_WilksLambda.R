#' @title Characteristic function of a linear combination of independent
#' LOG-TRANSFORMED WILK's LAMBDA distributed random variables
#'
#' @description
#' \code{cf_LogRV_WilksLambda(t, p, n, q, coef, niid)} evaluates characteristic function
#' of a linear combination (resp. convolution) of independent
#' LOG-TRANSFORMED WILK's LAMBDA distributed random variables.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param p vector of the dimension parameters \eqn{p = (p_1,...,p_N)}.
#' If empty, default value is \code{p = (1,...,1)}.
#' @param n vector of degrees of freedom of the Wishart matrices \eqn{E_i},
#' \eqn{n = (n_1,...,n_N)}. If empty, default value is \code{n = (1,...,1)}.
#' @param q vector of degrees of freedom of the Wishart matrices \eqn{H_i},
#' \eqn{q = (q_1,...,q_N)}. If empty, default value is \code{q = (1,...,1)}.
#' @param coef vector of the coefficients of the linear combination of the log-transformed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal.
#' If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient niid, such that \eqn{Z = Y + ... + Y}
#' is sum of niid random variables Y, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default \code{niid = 1}.
#'
#' @details
#' \code{cf_LogRV_WilksLambda} evaluates the characteristic function
#' of a random variable \eqn{Y  = coef_1*W_1 +...+ coef_N*W_N}, such that
#' \eqn{cf_Y(t) = cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t)}, where \eqn{cf_W_i(t)} is CF of \eqn{W_i = log(\Lambda_i)},
#' and each \eqn{\Lambda_i} has central WILK's LAMBDA distribution, \eqn{\Lambda_i ~ Lambda(p_i,n_i,q_i)},
#' for \eqn{i = 1,...,N}.
#' In particular \eqn{\Lambda_i = det(E_i)/det(E_i + H_i)}, with  \eqn{E_i} and \eqn{H _i}
#' being independent random matrices with central Wishart distributions,
#' \eqn{E_i ~  Wp_i(n_i,\Sigma_i)} and \eqn{H_i ~ Wp_i(q_i,\Sigma_i)} with \eqn{n_i >= p_i} and \eqn{\Sigma_i > 0}
#' (an unknown positive definite symmetric covariance matrix), for all \eqn{i = 1,...,N}.
#'
#' Each particular Lambda_i statistic can be used to test specific null hypothesis
#' (by measuring its effect expressed by the matrix \eqn{H_i}). In this case,
#' the null hypothesis is rejected for small values of the observed statistic \eqn{\Lambda_i},
#' or large value of \eqn{-log(\Lambda_i)}.
#'
#'  The central Wilks' distribution of \eqn{\Lambda_i ~ Lambda(p_i,n_i,q_i)},
#'  with \eqn{Lambda_i} in \eqn{(0,1)}, is defined by \eqn{\Lambda_i ~ Prod_{j=1}^p_i B_{i,j}},
#'  where \eqn{B_{i,j} ~ Beta{(n_i+1-j)/2, q_i/2)}}, i.e. \eqn{B_{i,j}} follow independent Beta distributions
#'  for all \eqn{i = 1,...,N} and \eqn{j = 1,...,p_i}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent
#' LOG-TRANSFORMED WILK's LAMBDA distributed random variables.
#'
#' @references
#' [1] Witkovsky, V., 2018. Exact distribution of selected multivariate test
#' criteria by numerical inversion of their characteristic functions.
#' arXiv preprint arXiv:1801.02248.
#'
#' @seealso For more details see WIKIPEDIA:
#' https://en.wikipedia.org/wiki/Wilks%27s_lambda_distribution.
#'
#' @family Continuous Probability Distribution
#'
#' @note Ver.: 04-Oct-2018 13:30:38 (consistent with Matlab CharFunTool v1.3.0, 19-Jul-2018 17:23:23).
#'
#' @example R/Examples/example_cf_LogRV_WilksLambda.R
#'
#' @export
#'
cf_LogRV_WilksLambda <- function(t, p, n, q, coef, niid) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(niid)) {
                niid <- vector()
        }
        if(missing(coef)) {
                coef <- vector()
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
        if(length(coef) == 0) {
                coef <- 1
        }
        if(length(niid) == 0) {
                niid <- 1
        }

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

        ## Characteristic function of a linear combination
        szt <- dim(t)
        t <- c(t)

        cf <- 1
        for(i in 1:length(coef)) {
                alpha <- (n[i] + 1 - (1:p[i])) / 2
                beta <- q[i] / 2
                cf <- cf * cf_LogRV_Beta(coef[i] * t, alpha, beta)
        }
        dim(cf) <- szt
        cf[t==0] <- 1

        if(length(niid) > 0) {
                if(length(niid) == 1) {
                        cf = cf ^ niid
                } else {
                        stop('niid should be a scalar (positive integer) value')
                }
        }

        return(cf)
}
