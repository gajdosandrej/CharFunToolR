#' @title The characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality  of mean vectors
#'
#' @description
#' Characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality of
#' mean vectors of q normal p-dimensional populations \eqn{(q>1)}. For more
#' details see Anderson (2003) and/or Marques and Coelho (2011).
#'
#' In particular, here we consider minus log-transformed LRT statistic
#' (Likelihood Ratio Test) or alternatively the minus log-transformed
#' modified LRT statistic.
#'
#'Let \eqn{X_k ~ N_p(mu_k,Sigma)} are p-dimensional random vectors with common
#'covariance matrix Sigma for all \eqn{k = 1,...,q}. We want to test the
#'hypothesis that the mean vectors mu_k are equal for all X_k, \eqn{k = 1,...,q}.
#'Then, the null hypothesis is given as
#'H0: \deqn{mu_1 = ... = mu_q},
#'i.e. the mean vectors are equal in all q populations. Here, the LRT test
#'statistic is given by
#'\deqn{LRT = Lambda^(n/2) = (det(E) / det(E+H))^(n/2)},
#'and the modified LRT is defined as
#'\deqn{MLRT = Lambda = det(E) / det(E+H)}
#' where \eqn{Lambda = det(E) / det(E+H)} with \eqn{E = sum_{k=1}^q sum_{j=1}^{n_k}
#' (X_{kj} - bar{X}_k)'*(X_{kj} - bar{X}_k)} with \eqn{E ~ Wishart(n-q,Sigma)},
#' and \deqn{H = sum_{k=1}^q (bar{X}_k - bar{X})'*(bar{X}_k - bar{X})} with H ~
#'  Wishart(q-1,Sigma) based on \eqn{n = n_1 + ... + n_q}  samples from the q
#' p-dimensional populations.
#'
#'
#' Under null hypothesis, distribution of the LRT statistic is
#' \eqn{LRT ~  prod_{j=1}^{q} (B_j)^{n/2}},
#' with \eqn{B_{j} ~  Beta((n-q-j+1)/2,(q-1)/2)}. Here we assume \eqn{n> p+q-1}.
#'
#'
#' @references
#' [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate
#' Statistical Analysis. New York: Wiley, 3rd Ed., 2003.
#'
#' [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
#' near-exact distribution theory for the most common likelihood ratio
#' test statistics used in Multivariate Analysis. Test, 2011, 20.1: 180-203.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param n total number of samples, \eqn{n = n_1 + ... + n_q } from the q p-dimensional populations.
#' If all sampes sizes are equal to n_1, then \eqn{n = n_1*q}  It is assumed that\eqn{n > p+q-1}.
#' @param p dimension of each of the q populations. It is assumed that each
#' population has the same (equal) dimension p.
#' @param q number of populations, \eqn{q>1}. If empty or missing, the default
#' value of q is \eqn{q = 2} (p is scalar).
#' @param type specifies the type of the LRT test statistic: \code{type = 'standard'}
#' with \eqn{W = -(n/2)*log(Lambda)}, or \code{type = 'modified'} with \eqn{W =
#' -log(Lambda)}. If type is empty or missing, default value is \code{type = 'modified'}.
#'
#' @return Characteristic function \code{cf} of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality of mean vectors.
#'
#' @note Ver.: 14-Sep-2022 15:02:57 (consistent with Matlab CharFunTool, 09-Dec-2018 09:23:48).
#'
#' @example R/Examples/example_cfTest_EqualityMeans.R
#'
#' @export


## ALGORITHM

cfTest_EqualityMeans<-function(t,n,p,q,type){

## CHECK THE INPUT PARAMETERS
if (missing(type)) {
        type <- vector()
        }
if (missing(q)){
        q <- vector()
        }

if (length(type)==0){
        type <- 'modified'
        }

if (length(q)==0){
        q <-2
                }

## Characteristic function of the null distribution
N <- p+q-1
if (n <= N){
        error('Sample size n is too small')
        }

ind   <- 1:p
alpha <- (n-q-ind+1)/2
beta  <- (q-1)/2

switch (type,
        "standard" =  cat(coef <- -n/2),
        "modified" = cat(coef <- -1),

        coef <- -1
)
cf   <- cf_LogRV_Beta(t,alpha,beta,coef)

return(cf)

}
