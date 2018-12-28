#' @title Characteristic function of a linear combination of \eqn{N} independent
#' LOG-TRANSFORMED MEANS-RATIO random variables
#'
#' @description
#' \code{cf_LogRV_MeansRatio(t, n, alpha, coef, niid)} evaluates characteristic function
#' of a linear combination (resp. convolution) of \eqn{N} independent LOG-TRANSFORMED MEANS-RATIO
#' random variables (RVs), \eqn{W_i = log(R_i)} for \eqn{i = 1,...,N} , where each \eqn{R_i = G_i/A_i} is a ratio
#' of \eqn{G_i} (the geometric mean) and \eqn{A_i} (the arithmetic mean) of a random sample
#' \eqn{X_{i,1},...,X_{i,n_i}} from a GAMMA distribution.
#'
#' That is, cf_LogRV_MeansRatio evaluates the characteristic function
#' of a random variable \eqn{Y  = coef_1*W_1 +...+ coef_N*W_N}, such that
#' \eqn{cf_Y(t) =  cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t)}, where \eqn{cf_W_i(t)}
#' is CF of \eqn{W_i = log(R_i)}, and  \eqn{R_i ~ MeansRatio(n_i,\alpha_i)},
#' where \eqn{MeansRatio(n_i,\alpha_i)} denotes the distribution of \eqn{R_i}, for \eqn{i = 1,...,N}.
#'
#' Here, we define \eqn{G_i = (X_{i,1} *...* X_{i,n_i})^(1/n_i)} and
#' \eqn{A_i = (X_{i,1} +...+ X_{i,n_i})/n_i}, with \eqn{X_{i,j} ~ Gamma(\alpha_i,\beta_i)}
#' for all \eqn{j = 1,...,n_i}, where \eqn{\alpha_i} is the shape parameter and \eqn{\beta_i}
#' is the rate parameter of the GAMMA distribution in the i-th random sample.
#'
#' Note that the \eqn{R_i} random variables are scale invariant,
#' so the distribution does not depend on the rate (scale) parameters \eqn{\beta_i}.
#'
#' The MeansRatio distribution of \eqn{R_i ~ MeansRatio(k_i,\alpha_i)},
#' with \eqn{R_i} in \eqn{(0,1)}, is defined by \eqn{R_i ~ (Prod_{j=1}^{k_i-1} B_{i,j})^(1/k_i)},
#' where \eqn{B_{i,j} ~ Beta{\alpha_i, j/k_i)}} for \eqn{j = 1,...,k_{i-1}},
#' i.e. \eqn{B_{i,j}} follow independent Beta distributions for all \eqn{i = 1,...,N} and
#' \eqn{j = 1,...,k_{i-1}}. Alternatively, \eqn{log(R_i) ~ (log(B_{i,1}) +...+ log(B_{i,k_{i-1}}))/k_i}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param n vector of sample size parameters \eqn{n = (n_1,...,n_N)}.
#' If empty, default value is \code{n = (1,...,1)}.
#' @param alpha vector of the 'shape' parameters \eqn{\alpha = (\alpha_1,...,\alpha_N)}.
#' If empty, default value is \code{alpha = (1,...,1)}.
#' @param coef vector of the coefficients of the linear combination of the
#' log-transformed random variables. If \code{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is
#' sum of \code{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N}
#' \eqn{coef(i) * log(X_i)} is independently and identically distributed
#' random variable. If empty, default value is \code{niid = 1}.
#'
#' @details
#'  The MeansRatio distribution is the distribution of the Bartlett's test
#'  statistitic for testing homogeneity of variances of \eqn{k} populations,
#'  based on equal sample sizes of size \eqn{m}, i.e. with \eqn{df = m-1}
#'  for all \eqn{j = 1,...,k}. Let \eqn{DF = k*df}. The Bartlett test statistic defined
#'  in Glaser (1976b) is
#'  \deqn{L = Prod((df*S^2_j/sigma^2)^(df/DF)) / (df/DF)*Sum(df*S^2_j/sigma^2) = Prod((S^2_j)^(1/k))/Sum(S^2_j)/k,}
#'  where \eqn{S^2_j = (1/(m-1))*Sum(X_jk - mean(X_jk))^2}.
#'  Notice that \eqn{df*S^2_j/sigma^2 ~ Chi^2_df = Gamma(df/2,1/2)} for all \eqn{j = 1,...,k}.
#'  The exact critical values can be calculated from the distribution of \eqn{W = -log(L)}
#'  by inverting the characteristic function \code{cf_LogRV_MeansRatio}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of N independent
#' LOG-TRANSFORMED MEANS-RATIO random variables.
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
#' https://en.wikipedia.org/wiki/Beta_distribution, \cr
#' https://en.wikipedia.org/wiki/Bartlett%27s_test.
#'
#' @family Continuous Probability Distribution
#'
#' @note Ver.: 05-Oct-2018 13:54:18 (consistent with Matlab CharFunTool v1.3.0, 17-Jun-2017 17:18:39).
#'
#' @example R/Examples/example_cf_LogRV_MeansRatio.R
#'
#' @export
#'
cf_LogRV_MeansRatio <- function(t, n, alpha, coef, niid) {
        ##CHECK THE INPUT PARAMETERS
        if(missing(niid)) {
                niid <- numeric()
        }
        if(missing(coef)) {
                coef <- vector()
        }
        if(missing(n)) {
                n <- vector()
        }
        if(missing(alpha)) {
                alpha <- vector()
        }

        if(length(alpha) == 0) {
                alpha <- 1
        }
        if(length(n) == 0) {
                n <- 1
        }
        if(length(coef) == 0) {
                coef <- 1
        }
        if(length(niid) == 0) {
                niid <- 1
        }

        ## Check size of the parameters
        l_max <- max(c(length(alpha), length(n), length(coef)))
        if (l_max > 1) {
                if (length(alpha) == 1) {
                        alpha <- rep(alpha, l_max)
                }
                if (length(n) == 1) {
                        n <- rep(n, l_max)
                }
                if (length(coef) == 1) {
                        coef <- rep(coef, l_max)
                }
                if ((any(lengths(list(
                        coef, alpha, n
                )) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        ## Characteristic function of a linear combination
        szt <- dim(t)
        t <- c(t)

        cf <- 1
        for(i in 1:length(coef)) {
                if(n[i] > 1) {
                        beta <- (1:(n[i]-1)) / n[i]
                        cf <- cf_LogRV_Beta(coef[i] * t, alpha[i], beta, 1 / n[i])
                }
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
