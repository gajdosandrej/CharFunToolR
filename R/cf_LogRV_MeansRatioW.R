#' @title Characteristic function of a linear combination of \eqn{N} independent
#' LOG-TRANSFORMED WEIGHTED MEANS-RATIO random variables
#'
#' @description
#' \code{cf_LogRV_MeansRatioW(t, n, alpha, weight, coef, niid)} evaluates characteristic function
#' of a linear combination (resp. convolution) of independent LOG-TRANSFORMED
#' WEIGHTED MEANS-RATIO random variables (RVs) \eqn{W_i = log(R_i)}, for \eqn{i = 1,...,N},
#' where each \eqn{R_i = G_i/A_i} is a ratio of the weighted geometric mean  \eqn{G_i}
#' and the (unweighted) arithmetic mean \eqn{A_i} of independent RVs \eqn{X_{i,1},...,X_{i,n_i}}
#' with GAMMA distributions, with the shape parameters \eqn{\alpha_{i,j}} and the (common) rate parameter
#' \eqn{\beta_i} for \eqn{i = 1,...,N} and \eqn{j = 1,...,n_i}.
#'
#' That is, \code{cf_LogRV_MeansRatioW} evaluates the characteristic function of a random
#' variable \eqn{Y = coef_1*W_1 +...+ coef_N*W_N}, such that \eqn{cf_Y(t) = cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t)},
#' where \eqn{cf_W_i(t)} is CF of \eqn{W_i = log(R_i)}, where \eqn{R_i} is the ratio statistic of the weighted geometric
#' mean and the arithmetic mean of independent gamma distributed RVs,
#' which distribution depends on the shape parameter \eqn{\alpha_{i,j}} and the vector
#' of weights \eqn{w_{i,j}}, for \eqn{i = 1,...,N} and \eqn{j = 1,...,n_i}.
#'
#' Here, for each fixed \eqn{i = 1,...,N}, the weighted geometric mean is defined
#' by \eqn{G_i = (X_{i,1}^w_{i,1} *...* X_{i,n_i}^w_{i,n_i})}, where the weights \eqn{w_{i,j}},
#' \eqn{j = 1,...,n_i}, are such that \eqn{w_{i,1} +...+ w_{i,n_i} = 1}, and the arithmetic mean is defined
#' by \eqn{A_i = (X_{i,1} +...+ X_{i,n_i})/n_i}.
#' The random variables  \eqn{X_{i,j}} are mutually independent with \eqn{X_{i,j} ~ \Gamma(\alpha_{i,j},\beta_i)}
#' for all \eqn{i = 1,...,N} and \eqn{j = 1,...,n_i}, where \eqn{\alpha_{i,j}} are the shape parameters and
#' \eqn{\beta_i} is the (common) rate parameter of the gamma distributions.
#'
#' Note that the weighted ratio random variables \eqn{R_i} are scale invariant,
#' so their distribution does not depend on the common rate (or scale) parameter \eqn{\beta_i}
#' for each \eqn{i = 1,...,N}.
#'
#' The distribution of the logarithm of the means ratio \eqn{log(R_i)} is defined
#' by its characteristic function, see e.g. Chao and Glaser (JASA 1978), which is
#' \deqn{cf_{log(R_i)}(t) = (n_i)^(1i*t) * ... \Gamma(sum(\alpha_{i,j}))/\Gamma(sum(\alpha_{i,j})+1i*t)) * ... Prod_{j=1}^n_i \Gamma(\alpha_{i,j}+1i*w_{i,j}*t)/\Gamma(\alpha_{i,j})},
#' for each \eqn{i = 1,...,N}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param n vector of sample size parameters \eqn{n = (n_1,...,n_N)}.
#' If empty, default value is \code{n = (1,...,1)}.
#' @param alpha list of weights vectors with shape parameters \eqn{alpha[[i]] = (\alpha_{i,1},...,alpha_{i,n_i})}.
#' If empty, default value is \code{alpha[[i]] = c(1,...,1)}.
#' @param weight list of weights vectors with \eqn{weights[[i]] = (w_{i,1},...,w_{i,n_i})}, for \eqn{i = 1,...,N}.
#' If empty, default value is \code{weights[[i]]  = c(1/k_i,...,1/k_i)}.
#' @param coef vector of the coefficients of the linear combination of the
#' log-transformed random variables. If \code{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is
#' sum of \code{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N}
#' \eqn{coef(i) * log(X_i)} is independently and identically distributed
#' random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of N independent
#' LOG-TRANSFORMED WEIGHTED MEANS-RATIO random variables.
#'
#' @references
#' Glaser, R. E. (1976a). The ratio of the geometric mean to the arithmetic
#' mean for a random sample from a gamma distribution.
#' \emph{Journal of the American Statistical Association}, 71(354), 480-487.
#'
#' Glaser, R. E. (1976b). Exact critical values for Bartlett's test
#' for homogeneity of variances. \emph{Journal of the American Statistical Association}, 71(354), 488-490.
#'
#' Chao, M. T., & Glaser, R. E. (1978). The exact distribution
#' of Bartlett's test statistic for homogeneity of variances with unequal sample sizes.
#' \emph{Journal of the American Statistical Association}, 73(362), 422-426.
#'
#' @seealso For more details see WIKIPEDIA:
#' https://en.wikipedia.org/wiki/Bartlett%27s_test.
#'
#' @family Continuous Probability Distribution
#'
#' @note Ver.: 05-Oct-2018 17:54:23 (consistent with Matlab CharFunTool v1.3.0, 17-Jun-2017 17:18:39).
#'
#' @example R/Examples/example_cf_LogRV_MeansRatioW.R
#'
#' @export
#'
cf_LogRV_MeansRatioW <- function(t, n, alpha, weight, coef, niid) {
        ##CHECK THE INPUT PARAMETERS
        if(missing(niid)) {
                niid <- numeric()
        }
        if(missing(coef)) {
                coef <- vector()
        }
        if(missing(weight)) {
                weight <- list()
        }
        if(missing(alpha)) {
                alpha <- list()
        }

        ## Check size of the parameters
        N <- length(n)

        if(length(coef) == 0){
                coef <- rep(1, N)
        }

        l_max <- max(c(length(n), length(coef)))
        if (l_max > 1) {
                if (length(n) == 1) {
                        n <- rep(n, l_max)
                }
                if (length(coef) == 1) {
                        coef <- rep(coef, l_max)
                }
                if ((any(lengths(list(
                        coef, n
                )) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        if(length(niid) == 0) {
                niid <- 1
        }

        if(length(alpha) == 0) {
                alpha <- list()
                for(i in 1:N) {
                        alpha[[i]] <- rep(1, n[i])
                }
        }

        if(length(weight) == 0) {
                weight <- list()
                for(i in 1:N) {
                        weight[[i]] <- rep(1, n[i]) / n[i]
                }
        }

        if(!is.list(alpha)) {
                alpha <- as.list(alpha)
        }

        if(!is.list(weight)) {
                weight <- as.list(weight)
        }

        for(i in 1:N) {
                aux <- rep(1, n[i])
                l_max <- max(c(length(alpha[[i]]), length(weight[[i]]), length(aux)))
                if (l_max > 1) {
                        if (length(alpha[[i]]) == 1) {
                                alpha[[i]] <- rep(alpha[[i]], l_max)
                        }
                        if (length(weight[[i]]) == 1) {
                                weight[[i]] <- rep(weight[[i]], l_max)
                        }
                        if (length(aux) == 1) {
                                aux <- rep(aux, l_max)
                        }
                        if ((any(lengths(list(
                                aux, weight[[i]], alpha[[i]]
                        )) < l_max))) {
                                stop("Input size mismatch.")
                        }
                }

        }

        ## Characteristic function of a linear combination
        szt <- dim(t)
        t <- c(t)

        cf <- 1
        for(i in 1:N) {
                cf <- cf * cf_LogR(coef[[i]] * t, n[i], alpha[[i]], weight[[i]])
        }

        dim(cf) <- szt
        cf[t==0] <- 1

        if(length(niid) > 0) {
                if(length(niid) == 1) {
                        cf <- cf ^ niid
                } else {
                        stop('niid should be a scalar (positive integer) value.')
                }
        }

        return(cf)
}



## Function cf_LogR
cf_LogR <- function(t, k, alpha, weight) {
        # Auxiliary function. Calculates the characteristic function of the random
        # variable log(R), where R is a ratio of weighted geometric mean and the
        # unweighted arithmetic mean of independent gamma random variables with
        # their shape parameters given by the vector alpha.

        A <- sum(unlist(alpha))
        cf <- (1i * t) * log(k) + GammaLog(A) - GammaLog(A +1i * t)
        for(j in 1:k) {
                cf <- cf + GammaLog(alpha[[j]] + 1i * t * weight[[j]]) - GammaLog(alpha[[j]])
        }

        cf <- exp(cf)
}

