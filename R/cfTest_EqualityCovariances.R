#' @title The characteristic function of the the exact null distribution of the multivariate analysis (MVA) test statistic.
#'
#' @description
#' Characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality of
#' covariance matrices of q normal p-dimensional populations \eqn{(q>1)}. For
#' more details see Anderson (2003) and/or Marques and Coelho (2011).
#'
#' In particular, here we consider minus log-transformed LRT statistic
#' (Likelihood Ratio Test) or alternatively the minus log-transformed
#' modified LRT statistic.
#'
#'Let \eqn{X_k ~ N_p(mu_k,Sigma_k)} are p-dimensional random vectors, for \eqn{k =
#'1,...,q}. We want to test the hypothesis that the covariance matrix Sigma
#'is common for all X_k, \eqn{k = 1,...,q}. Then, the null hypothesis is given as
#'H0: \deqn{Sigma_1 = ...sigma_q}
#'i.e. the covariance matrices are equal in all q populations. Here, the
#'LRT test statistic is given by
#'\deqn{LRT = Lambda^(n/2) = (q^(p*q) * prod(det(S_k)) / (det(S))^q)^(n/2)},
#' and the modified LRT is defined as
#' \eqn{MLRT = Lambda = q^(p*q) * prod(det(S_k)) / (det(S))^q},
#' where \eqn{Lambda = q^(p*q) * prod(det(S_k)) / (det(S))^q} where S_k are MLEs
#' of Sigma_k, for \eqn{k = 1,...,q}, and \eqn{S = S_1 + ... + S_q}, based on n samples
#' from each of the the q p-dimensional populations.
#'
#' Under null hypothesis, distribution of the test LRT statistic is
#' \eqn{LRT ~  prod_{k=1}^{q} prod_{j=1}^{p} (B_{jk})^{n/2}},
#' with \eqn{B_{jk} ~ Beta((n-j)/2,(j*(q-1)+2*k-1-q)/2)}, and we set \eqn{B_{11} = 1}
#' for \eqn{j=k=1}. Here we assume that \eqn{n > p}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic Lambda, say \eqn{W = -log(LRT)} is
#' given by \code{cf = function(t) cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q))},
#' where \eqn{k = [1*o,...,q*o]} with p-dimensional vector of ones \eqn{o = [1,...,1]}
#' and \eqn{j = [j_1,...,j_q]} with \eqn{j_k = 1:p}.
#' Similarly, the exact characteristic function of the null distribution of
#' minus log-transformed modified LRT statistic, say \eqn{W = -log(MLRT)} is
#' \code{cf = function(t) cf_LogRV_Beta(-t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q))},
#' where \eqn{k = [1*o,...,q*o]} with p-dimensional vector of ones \eqn{o = [1,...,1]}
#' and \eqn{j = [j_1,...,j_q]} with \eqn{j_k = 1:p}.
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
#' @param n number of samples from each of the q p-dimensional populations.
#' It is assumed that\eqn{n > p}.
#' @param p dimension of each of the q populations. It is assumed that each
#' population has the same (equal) dimension p.
#' @param q number of populations, \eqn{q>1}. If empty or missing, the default
#' value of q is \eqn{q = 2} (p is scalar).
#' @param type specifies the type of the LRT test statistic: \code{type = 'standard'}
#' with \eqn{W = -(n/2)*log(Lambda)}, or \code{type = 'modified'} with \eqn{W =
#' -log(Lambda)}. If type is empty or missing, default value is \code{type = 'modified'}.
#'
#' @return Characteristic function \code{cf} of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality of covariance matrices
#'
#' @note Ver.: 14-Sep-2022 14:26:55 (consistent with Matlab CharFunTool, 09-Dec-2018 09:23:48).
#'
#' @example R/Examples/example_cfTest_EqualityCovariances.R
#'
#' @export
#'
#  ALGORITHM
cfTest_EqualityCovariances<-function(t,n,p,q,type){

## CHECK THE INPUT PARAMETERS


if (missing(type)){
        type <- vector}
if (missing(q)){
        q <- vector()}

if (length(type)==0){
type = 'modified'
}

if (length(q)==0){
q <- 2
}

## Characteristic function of the null distribution
if (n <= p){
error('Sample size n is too small')
}

alpha <- rep(0,p*q)
beta  <- rep(0,p*q)
ind   <- 0
for (k in 1:q){
for (j  in 1:p){
ind <- ind +1
alpha[ind] <- (n - j) / 2
beta[ind]  <- (j*(q-1) + 2*k - 1 - q) / (2*q)
}
}

# For k=j=1 the coefficient beta=0, hence log(Beta(alpha,beta))=0.
alpha <- alpha[2:length(alpha)]
beta  <- beta[2:length(beta)]

switch (type,
        "standard" =  cat(coef <- -n/2),
        "modified" = cat(coef <- -1),

        coef <- -1
)


cf   <- cf_LogRV_Beta(t,alpha,beta,coef)

return(cf)
}

