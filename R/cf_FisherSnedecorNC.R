#' @title Characteristic function of a linear combination
#'  of independent non-central Fisher-Snedecor random variables
#'
#' @description
#' \code{cf_FisherSnedecorNC(t, df1, df2, delta, coef, niid, tol)} evaluates characteristic function
#' of a linear combination (resp. convolution) of independent non-central Fisher-Snedecor random variables,
#' with distributions \eqn{F(df1_i,df2_i,\delta_i)}.
#'
#' That is, cf_FisherSnedecorNC evaluates the characteristic function
#' \eqn{cf(t)} of  \eqn{Y = coef_i*X_1 +...+ coef_N*X_N}, where \eqn{X_i ~ F(df1_i,df2_i,\delta_i)}
#' are inedependent RVs, with \eqn{df1_i} and \eqn{df2_i} degrees of freedom,
#' and the noncentrality parameters \eqn{\delta_i >0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of X ~ F(df1,df2,delta) is Poisson mixture of the CFs of the scaled central F RVs
#' of the form
#' \deqn{cf(t) = cf_FisherSnedecorNC(t,df1,df2,\delta) = exp(-\delta/2) sum_{j=1}^Inf (\delta/2)^j/j! * cf_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2),}
#' where \eqn{cf_FisherSnedecor(t,df1,df2)} are the CFs of central F RVs with parameters \eqn{df1} and \eqn{df2}.
#' Hence,the characteristic function of \eqn{Y = coef(1)*Y1 + ... + coef(N)*YN is cf_Y(t) = cf_Y1(coef(1)*t) * ... *cf_YN(coef(N)*t)},
#' where \eqn{cf_Yi(t)} is evaluated with the parameters \eqn{df1_i}, \eqn{df2_i}, and \eqn{\delta_i}.
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
#' of independent non-central Fisher-Snedecor random variables.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Noncentral_F-distribution}.
#'
#' @family Continuous Probability Distribution
#'
#' @example R/Examples/example_cf_FisherSnedecorNC.R
#'
#' @export
#'
cf_FisherSnedecorNC <- function(t, df1, df2, delta, coef, niid, tol) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(df1)) {
    df1 <- vector()
  }
  if(missing(df2)) {
    df2 <- vector()
  }
  if(missing(delta)) {
    delta <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- numeric()
  }
  if(missing(tol)) {
    tol <- numeric()
  }

  if(length(df1) == 0){
    df1 <- 1
  }
  if(length(df2) == 0){
    df2 <- 1
  }
  if(length(delta) == 0){
    delta <- 0
  }
  if(length(coef) == 0){
    coef <- 1
  }
  if(length(tol) == 0){
    tol <- 1e-12
  }

  ## SET THE COMMON SIZE of the parameters
  l_max <- max(c(length(df1), length(df2), length(delta), length(coef)))
  if (l_max > 1) {
    if (length(df1) == 1) {
      df1 <- rep(df1, l_max)
    }
    if (length(df2) == 1) {
      df2 <- rep(df2, l_max)
    }
    if (length(delta) == 1) {
      delta <- rep(delta, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(
      coef, df1, df2, delta
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function of a linear combination of independent nc F RVs
  szt <- dim(t)
  t <- c(t)
  cf <- rep(1, length(t))
  for(i in 1:length(coef)) {
    cf <- cf * cf_ncFisherSnedecor(coef[i] * t, df1[i], df2[i], delta[i], tol)
  }
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

## Function funCF
cf_ncFisherSnedecor <- function(t, df1, df2, delta, tol) {
  # cf_ncFisherSnedecor Characteristic function of the distribution of
  # the of the non-central Fisher-Snedecor distribution with df1 and df2
  # degrees of freedom and the non-centrality parameter delta > 0.
  f <- 0
  delta <- delta / 2
  if(delta == 0){
    # Deal with the central distribution
    f <- cf_FisherSnedecor(t, df1, df2)
  } else if(delta > 0) {
    # Sum the Poisson series of CFs of independent central F RVs,
    # poisspdf(j,delta).*cf_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2)
    j0 <- floor(delta / 2)
    p0 <- exp(-delta + j0 * log(delta) - log(gamma(j0 + 1)))
    f <- f + p0 * cf_FisherSnedecor(t * (df1 + 2 * j0) / df1, df1 + 2 * j0, df2)
    p <- p0
    j <- j0 - 1
    while(j >= 0 && p > tol) {
      p <- p *  (j + 1) / delta
      f <- f + p * cf_FisherSnedecor(t* (df1 + 2 * j) / df1, df1 + 2 * j, df2)
      j <- j - 1
    }
    p <- p0
    j <- j0 + 1
    i <- 0
    while(p > tol && i <= 5000) {
      p <- p * delta / j
      f <- f + p * cf_FisherSnedecor(t * (df1 + 2 * j) / df1, df1 + 2 * j, df2)
      j <- j + 1
      i <- i + 1
    }
    if(i == 5000) {
      warning("No convergence")
    }
    return(f)
  } else {
    stop("delta should be nonnegative")
  }
}
