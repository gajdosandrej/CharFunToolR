#' @title The characteristic function of the the exact null distribution of the BARTLETT's test statistic

#' @description
#' \code{cfTest_Bartlett(t,df)} evaluates the characteristic function of the exact null distribution of the BARTLETT's
#' test statistic for testing hypothesis about homogeneity of variances
#' of k normal populations, specified by the vector of the degrees of
#' freedom \deqn{df = [df_1,...,df_k]}.
#'
#' The characteristic function of the the exact null distribution of the
#' BARTLETT's test statistic is given by
#' \deqn{cf(t) = exp(1i*t*C/B).*
#' k.^(1i*t*D/B).*(gamma(D/2)/gamma(D/2-1i*t*D/B)).*
#' prod(gamma(df/2-1i*t*df/B)./gamma(df/2))} where \eqn{B = 1+1/(3*(k-1))*(sum(1./df)-1/D)} and \eqn{C = D*log(k*prod(W.^W))};
#' with \eqn{W = df/D}, and \eqn{D = sum(df)}.
#'
#'
#'
#' @references
#' [1]  Glaser, R. E. (1976a). The ratio of the geometric mean to the arithmetic
#' mean for a random sample from a gamma distribution. Journal of the
#' American Statistical Association, 71(354), 480-487.
#'
#' [2] Glaser, R. E. (1976b). Exact critical values for Bartlett's test for
#' homogeneity of variances. Journal of the American Statistical
#' Association, 71(354), 488-490.
#'
#' [3] Chao, M. T., & Glaser, R. E. (1978). The exact distribution of
#' Bartlett's test statistic for homogeneity of variances with unequal
#' sample sizes. Journal of the American Statistical Association, 73(362), 422-426.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Bartlett%27s_test},
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param df k-dimensional vector of degrees of freedom \eqn{df = [df_1,...,df_k]}.
#'
#' @return Characteristic function \eqn{cf(t)} of the  exact null distribution of the BARTLETT's test statistic.
#'
#' @note Ver.: 14-Sep-2022 11:42:43 (consistent with Matlab CharFunTool, 09-Dec-2018 09:42:32).
#'
#' @example R/Examples/example_cfTest_Bartlett.R
#'
#' @export
#'
## Algorithm
cfTest_Bartlett<-function(t,df){

szt <- dim(t)
t   <- t(t)
df  <- t(df);
k   <- length(df)
D   <- sum(df)
W   <- df/D;
C   <- D*log(k*prod(W^W))
B   <- 1 + (1/(3*(k-1))) * (sum(1/df)-1/D)

cf  <- GammaLog(D/2)
cf  <- cf + (1i*t*C/B) - log(k)*(1i*t*D/B) - GammaLog(D/2-1i*t*D/B)
for ( j in 1:k) {
        cf <- cf + GammaLog(df[j]/2-1i*t*df[j]/B) - GammaLog(df[j]/2)

}
cf<-exp(cf)
cf[t==0] <- 1
dim(cf)  <-szt

return(cf)
}
