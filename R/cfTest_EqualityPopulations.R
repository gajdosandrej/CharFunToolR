#' @title The characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality of q normal p-dimensional populations.
#'
#' @description
#' Characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing equality of q
#' normal p-dimensional populations. In fact, we simultaneously test
#' equality of means and equality of covariance matrices of q normal
#' p-dimensional populations. For more details see Anderson (2003) and/or
#' Marques and Coelho (2011).
#'
#' In particular, here we consider minus log-transformed LRT statistic
#' (Likelihood Ratio Test) or alternatively the minus log-transformed
#' modified LRT statistic.
#'
#'Let \eqn{X_k ~ N_p(mu_k,Sigma_k)}  for \eqn{k = 1,...,q}. We want to test the
#'hypothesis that the q normal populations are equally distributed. That
#'is, we want to test that the mean vectors mu_k are equal for all \eqn{k=
#'1,...,q}, as well as the covariance matrices Sigma_k are equal for all \eqn{k
#'= 1,...,q}. Then, the null hypothesis is given as
#'H0: \deqn{mu_1 = ... = mu_q} and \eqn{Sigma_1 = ... = Sigma_k}.
#' Here, the null hypothesis H0 and the LRT statistic can be decomposed:
#' \eqn{LRT = LRT_Means * LRT_Covariances}
#'where (first) LRT_Covariances represents the LRT for testing equality of
#'covariance matrices of given q normal populations, and (second)
#'LRT_Means represents (conditionally) the LRT for testing equality of
#'means of given q normal populations.
#'
#' Under null hypothesis, distribution of LRT_Covariances and LRT_Means
#' are independent, and the distribution of the compound LRT statistic is
#' \deqn{LRT =  Lambda^(n/2) = LRT_Means * LRT_Covariances = (prod_{k=1}^q prod_{j=1}^{p} (B_{jk})^(n/2)) * (prod_{j=1}^{p} (B_j)^(q*n/2))},
#' and the modified LRT is defined as
#'\deqn{ MLRT = Lambda = MLRT_Means * MLRT_Covariances = (prod_{k=1}^q prod_{j=1}^{p} B_{jk}) * (prod_{j=1}^{p} (B_j)^q)}
#'
#'  where \deqn{Lambda = (prod_{k=1}^q prod_{j=1}^p B_{jk})*(prod_{j=1}^{p} B_j^q)}
#'  where  B_{jk} and B_j are mutually independent beta distributed
#'  random variables. Here we assume that n is equal sample size for each
#'  sample, \eqn{k = 1,...,q}, \eqn{n > p}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic Lambda, say W = -log(LRT) is
#' given by
#'  \code{cf = function(t) cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q)) .* cf_LogRV_Beta(-(n*q/2)*t, ((n-1)*q-i+1)/2, (q-1)/2)},
#'   where \eqn{i = [1, 2, ..., p]'}, \eqn{k = [1*o,...,q*o]} with p-dimensional vector
#'   of ones \eqn{o = [1,...,1]}  and \eqn{j = [j_1,...,j_q]} with \eqn{j_k = 1:p}.
#'   Similarly, the exact characteristic function of the null distribution of
#'   minus log-transformed modified LRT statistic, say \eqn{W = -log(MLRT)} is
#'     \code{cf = function(t) cf_LogRV_Beta(-t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q)) .* cf_LogRV_Beta(-q*t, ((n-1)*q-i+1)/2, (q-1)/2)},
#'    where \eqn{i = [1, 2, ..., p]'}, \eqn{k = [1*o,...,q*o]} with p-dimensional vector
#'    of ones \eqn{o = [1,...,1]}  and \eqn{j = [j_1,...,j_q]} with \eqn{j_k = 1:p}.
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
#' @param n number of samples, from each of the q p-dimensional populations.
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
#' multivariate analysis (MVA) test statistic for testing equality of q normal p-dimensional populations.
#'
#' @note Ver.: 14-Sep-2022 15:32:11 (consistent with Matlab CharFunTool, 09-Dec-2018 09:23:48).
#'
#' @example R/Examples/example_cfTest_EqualityPopulations.R
#'
#' @export


## ALGORITHM
cfTest_EqualityPopulations<-function(t,n,p,q,type){

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

        if (n <= p){
                error('Sample size n is too small')
        }

        ind   <- 1:p
        alpha <- rep(0,p*(q+1))
        beta  <- rep(0,p*(q+1))
        coef<-rep(1,p*(q+1))

        alpha[1:p] <- ((n-1)*q-ind+1)/2
        beta[1:p]  <- (q-1)/2
        coef[1:p]  <- q*coef[1:p]

        ind    <- p
        for (k  in 1:q){
        for (j in 1:p){
        ind <- ind +1
        alpha[ind] <- (n - j) / 2
        beta[ind]  <- (j*(q-1) + 2*k - 1 - q) / (2*q)
        }
        }


        switch (type,
                "standard" =  cat(coef <- -(n/2)*coef),
                "modified" = cat(coef <- -coef),

                coef <- -coef
        )
        cf   <- cf_LogRV_Beta(t,alpha,beta,coef)

        return(cf)

}


