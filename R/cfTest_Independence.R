#' @title The characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing independence q normal populations.
#'
#' @description
#' Characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing independence q
#' normal populations. For more details see Anderson (2003) and/or Marques
#' and Coelho (2011).
#'
#' In particular, here we consider minus log-transformed LRT statistic
#' (Likelihood Ratio Test) or alternatively the minus log-transformed
#' modified LRT statistic.
#'
#'Let \eqn{X_k ~ N_p(mu_k,Sigma_k)}  are p_k dimensional random vectors, for \eqn{k = 1,...,q}. Let us denote X = [X_1,...,X_q] and assume X ~N_p(mu,Sigma). Then, the null hypothesis is given as
#'H0: \deqn{Sigma = diag(Sigma_1,...,Sigma_q)},
#'i.e. the off-diagonal blocks of Sigma are blocks of zeros. Here, the LRT
#'test statistic is given by
#'\deqn{LRT = Lambda^(n/2) = (det(S)/prod(det(S_k)))^(n/2)},
#'and the modified LRT is defined as
#'\deqn{MLRT = Lambda = det(S)/prod(det(S_k))},
#'where \eqn{Lambda = det(S)/prod(det(S_k))} with S is MLE of Sigma, and S_k are
#'MLEs of Sigma_k, for \eqn{k = 1,...,q}, based on n samples from the compound
#'vector \eqn{X = [X_1,...,X_q]}. Here we assume that \eqn{n > p = sum(p_k)}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic Lambda, say\eqn{W = -log(LRT)} is
#' given by
#'  \code{cf = function(t) cf_LogRV_Beta(-(n/2)*t, (n-qk-j)/2, qk/2)},
#'  where \eqn{qk = [qk_1,...,qk_q]} with \deqn{qk_k = p_{k+1} + ... + p_m}, and \eqn{j = [j_1,...,j_m]} with \eqn{j_k = 1:p_k}.
#' Similarly, the exact characteristic function of the null distribution of
#' minus log-transformed modified LRT statistic, say \eqn{W = -log(MLRT)} is
#' \code{cf = function(t) cf_LogRV_Beta(-t, (n-qk-j)/2, qk/2)},
#' where \eqn{qk = [qk_1,...,qk_m]} with \eqn{qk_k = p_{k+1} + ... + p_m}, and \eqn{j =[j_1,...,j_m]} with \eqn{j_k = 1:p_k}.
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
#' @param n number of samples, from the compound vector \eqn{X = [X_1,...,X_q]}. It is assumed that \eqn{n > length(X)}.
#' @param p q-dimensional vector of dimensions of the vectors \eqn{X_1,...,X_q}.
#' If p is scalar value, then the algorithm assumes equal dimension
#' p for all q vectors \eqn{X_1,...,X_q}.
#' @param q number of populations, \eqn{q>1}. If p is vector with dimension\eqn{ >= 2},
#' then q must be equal to the value \eqn{q = length(p)}. If empty or
#' missing, the default value of q is \eqn{q = length(p)} (p is vector)
#' or \eqn{q = 2} (p is scalar).
#' @param type specifies the type of the LRT test statistic: \code{type = 'standard'}
#' with \eqn{W = -(n/2)*log(Lambda)}, or \code{type = 'modified'} with \eqn{W =
#' -log(Lambda)}. If type is empty or missing, default value is \code{type = 'modified'}.
#'
#' @return Characteristic function \code{cf} of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing independence q normal populations.
#'
#' @note Ver.: 14-Sep-2022 15:56:09 (consistent with Matlab CharFunTool, 09-Dec-2018 09:23:48).
#'
#' @example R/Examples/example_cfTest_EqualityIndepedence.R
#'
#' @export


## ALGORITHM
cfTest_Indepedence<-function(t,n,p,q,type){

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
                if(length(p)>1){
                q <-length(p)
                }
                else{
                        q<-2
                        p<-c(p,p)
                }
        }
        else{
                if(length(p)>1&&length(p)!=q){
                 warning('Dimensions mismatch. Dimensions of p must be equal to q')
                        q<-length(p)
                }
                else if(length(p)==1){
                        p<-p*rep(1,q)
                }
        }

        ## Characteristic function of the null distribution

        csp<-cumsum(p)
        if (n <= csp[q]){
                error('Sample size n is too small')
        }


        alpha <- rep(0,csp[q-1])
        beta  <- rep(0,csp[q-1])
        ind   <- 0


        for (k  in 1:(q-1)){
                qk<-csp[q]-csp[k]
                for (j in 1:p[k]){
                        ind <- ind +1
                        alpha[ind] <- (n - qk-j) / 2
                        beta[ind]  <- qk/2
                }
        }


        switch (type,
                "standard" =  cat(coef <- -n/2),
                "modified" = cat(coef <- -1),

                coef <- -1
        )
        cf   <- cf_LogRV_Beta(t,alpha,beta,coef)

        return(cf)

}
