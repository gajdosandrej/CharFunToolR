#' @title Characteristic function of a linear combination
#' of independent  LOG-TRANSFORMED FISHER-SNEDECOR F random variables
#'
#' @description
#' \code{cf_LogRV_FisherSnedecor(t, df1, df2, coef, niid)} evaluates characteristic function
#' of a linear combination (resp. convolution) of independent  LOG-TRANSFORMED FISHER-SNEDECOR F
#' random variables (RVs) \eqn{log(X)}, where \eqn{X ~ F(df1,df2)} has the FISHER-SNEDECOR F distribution
#' with \eqn{df1 > 0} and \eqn{df2 > 0} degrees of freedom.
#'
#' That is, \code{cf_LogRV_FisherSnedecor} evaluates the characteristic function
#' \eqn{cf(t) of  Y = coef_1*log(X_1) +...+ coef_N*log(X_N)}, where X_i ~ F(df1_i,df2_i),
#' with degrees of freedom \eqn{df1_i} and \eqn{df2_i}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y = log(X)}, with \eqn{X ~ F(df1,df2,\lambda)},
#' where \eqn{\lambda} is the non-centrality parameter, is defined by
#' \eqn{cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t))}.
#' That is, the characteristic function can be derived from expression
#' for the r-th moment of \eqn{X}, \eqn{E(X^r)} by using \eqn{(1i*t)} instead of \eqn{r}.
#' In particular, the characteristic function of \eqn{Y = log(X)} is defined by
#' \deqn{cf_Y(t) = (df2/df1)^(1i*t) * gamma(df1/2 + 1i*t) / gamma(df1/2) * gamma(df2/2 - 1i*t) / gamma(df2/2) * 1F1(-1i*t;df1/2;-\lambda/2),}
#' where \eqn{1F1(a;b;z)} is the hypergeometric function confluent hypergeometric function,
#' also known as the Kummer function \eqn{M(a,b,z)}.
#'
#' Hence,the characteristic function of \eqn{Y  = coef(1)*Y1 + ... + coef(N)*YN}
#' is  \deqn{cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t),} where \eqn{cf_Yi(t)}
#' is evaluated with the parameters \eqn{df1(i)}, \eqn{df2(i)}, and \eqn{\lambda(i)}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param df1 vector of the  degrees of freedom \code{df1 > 0}. If empty, default value is \code{df1 = 1}.
#' @param df2 vector of the  degrees of freedom \code{df2 > 0}. If empty, default value is \code{df2 = 1}.
#' @param coef vector of the coefficients of the linear combination of the Beta distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#' @param tol tolerance factor for selecting the Poisson weights, i.e. such that \eqn{PoissProb > tol}.
#' If empty, default value is \code{tol = 1e-12}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent  LOG-TRANSFORMED FISHER-SNEDECOR F random variables.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/F-distribution}.
#'
#' @family Continuous Probability Distribution
#'
#' @example R/Examples/example_cf_LogRV_FisherSnedecor.R
#'
#' @export
#'
cf_LogRV_FisherSnedecor <- function(t, df1, df2, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(df1)) {
    df1 <- vector()
  }
  if(missing(df2)) {
    df2 <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- numeric()
  }

  if(length(df1) == 0) {
    df1 <- 1
  }
  if(length(df2) == 0) {
    df2 <- 1
  }
  if(length(coef) == 0) {
    coef <- 1
  }
  if(length(niid) == 0) {
    niid <- 1
  }

  ## Check size of the parameters
  l_max <- max(c(length(df1), length(df2), length(coef)))
  if (l_max > 1) {
    if (length(df1) == 1) {
      df1 <- rep(df1, l_max)
    }
    if (length(df2) == 1) {
      df2 <- rep(df2, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(
      coef, df1, df2
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function of a linear combination
  szt <- dim(t)
  t <- c(t)
  aux <- 1i * t %*% t(coef)
  aux <- GammaLog(t(t(aux) + (df1 / 2)))  + GammaLog(t(-t(aux) + (df2 / 2))) - rep(1, length(t)) %*% t(GammaLog(df1 / 2) + GammaLog(df2 / 2)) + t(t(aux) * (log(df2 / df1)))
  cf <- apply(exp(aux), 1, prod)

  dim(cf) <- szt
  cf[t==0] <- 1

  if(length(niid) > 0) {
    if(length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value")
    }
  }

  return(cf)
}
