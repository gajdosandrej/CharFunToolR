#' @title The characteristic function of the the exact null distribution of the multivariate analysis (MVA) test statistic
#'
#' @description
#' Characteristic function of the exact null distribution of the
#' multivariate analysis (MVA) test statistic for testing the shape of the
#' unknown covariance matrix Sigma of p-dimensional normal distribution.
#' Here it is assumed that under null hypothesis the covariance matric has
#' the compound symmetry structure, i.e. \eqn{Sigma = sigma^2*((1-rho)*I+rho*J)},
#' where I is the p-dimensional identity matrix and J is pxp matrix of
#' ones. For more details see Anderson (2003), Marques and Coelho (2011),
#' Nagar and Castaneda(2006), Marques and Coelho (2013 and 2018).
#'
#' In particular, here we consider minus log-transformed LRT statistic
#' (Likelihood Ratio Test) or alternatively the minus log-transformed
#' modified LRT statistic.
#'
#'Let \eqn{X ~ N_p(mu,Sigma)} with \eqn{p >= 2}. Then the null hypothesis is given as
#'H0: \deqn{Sigma = sigma^2*((1-rho)*I_p + rho*J_p) (sigma^2>0} unspecified).
#'Here, the LRT test statistic is given by
#'\deqn{LRT = Lambda^(n/2) = ((p-1)^(p-1) * p^p * det(S) / ...
#'                (trace(J*S) * trace((p*I-J)*S)^(p-1))^(n/2)},
#'and the modified LRT is defined as
#' \deqn{MLRT = Lambda = (p-1)^(p-1) * p^p * det(S) / ...
#'                (trace(J*S) * trace((p*I-J)*S)^(p-1)},
#'where S is MLE of Sigma based on sample size n from \eqn{X ~ N_p(mu,Sigma)}.
#'
#' Under null hypothesis, distribution of the test statistic LRT is
#' \eqn{LRT ~  prod_{j=2}^{p} (B_j)^{n/2}},
#' with \eqn{B_j ~ Beta((n-j+1)/2,(j-2)/(p-1) + (j-1)/2)}, where \eqn{j = [2,...,p]}.
#' Here we assume that \eqn{n > p}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic Lambda, say \eqn{W = -log(LRT)} is
#' given by \code{cf = function(t) cf_LogRV_Beta(-(n/2)*t, (n-j+1)/2, (j-2)/(p-1) + (j-1)/2)},
#' where \eqn{j = [2,..., p]'}.
#' Similarly, the exact characteristic function of the null distribution of
#' minus log-transformed modified LRT statistic, say \eqn{W = -log(MLRT)} is
#' \code{cf = function(t) cf_LogRV_Beta(-t, (n-j+1)/2, (j-2)/(p-1) + (j-1)/2)},
#' where \eqn{j = [2,..., p]'}.
#'
#' @references
#' [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate
#' Statistical Analysis. New York: Wiley, 3rd Ed., 2003.
#'
#' [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
#' near-exact distribution theory for the most common likelihood ratio
#' test statistics used in Multivariate Analysis. Test, 2011, 20.1: 180-203.
#'
#' [3] NAGAR, D.K., CASTANEDA, E.M. Distribution and percentage
#' points of LRC for testing multisample compound symmetry in the
#' bivariate and the trivariate cases. Metron 64, no. 2 (2006): 217-238.
#'
#' [4] MARQUES, Filipe J.; COELHO, Carlos A. The multisample block-diagonal
#' equicorrelation and equivariance test, AIP Conference Proceedings 1558 (2013) pp. 793-796.
#'
#' [5] MARQUES, Filipe J.; COELHO, Carlos A. Testing simultaneously
#' different covariance block diagonal structures - the multisample
#' case. 2018 Working paper.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param n the sample size from the p-dimensional population. It is assumed
#' that \eqn{n > p}.
#' @param p dimension of the population \code{(p>=2)}.
#' @param type specifies the type of the LRT test statistic: \code{type = 'standard'}
#' with \eqn{W = -(n/2)*log(Lambda)}, or \code{type = 'modified'} with \eqn{W =
#' -log(Lambda)}. If type is empty or missing, default value is \code{type = 'modified'}.
#'
#' @return Characteristic function \eqn{cf(t)} of the exact null distribution of the multivariate analysis (MVA) test statistic.
#'
#' @note Ver.: 14-Sep-2022 13:59:03 (consistent with Matlab CharFunTool, 09-Dec-2018 09:23:48).
#'
#' @example R/Examples/example_cfTest_CompundSymmetry.R
#'
#' @export
#'
##ALGORITHM
cfTest_CompoundSymmetry<-function(t,n,p,type){

## CHECK THE INPUT PARAMETERS

if (missing(type)){
        type <- vector()
}

if (length(type)==0){
type <- 'modified'
}

## Characteristic function of the null distribution
if (n <= p ){
error('Sample size n is too small')
}

ind   <- 2:p
alpha <- (n-ind+1)/2
beta  <- (ind-2)/(p-1) + (ind-1)/2

switch (type,
    "standard" =  cat(coef <- -n/2),
    "modified" = cat(coef <- -1),

        coef <- -1
)

cf   <-cf_LogRV_Beta(t,alpha,beta,coef)

return(cf)

}

