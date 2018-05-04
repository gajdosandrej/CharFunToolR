#' @title Characteristic function of a linear combination of independent LOG-TRANSFORMED non-central ChiSquare random variables
#'
#' @description
#' \code{cf_LogRV_ChiSquareNC(t, df, delta, coef, niid, tol)} evaluates characteristic function
#' of a linear combination (resp. convolution) of independent LOG-TRANSFORMED non-central ChiSquare random variables,
#' with distributions \eqn{ChiSquare(df_i,\delta_i)}.
#'
#' That is, \code{cf_LogRV_ChiSquareNC} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = coef_i*log(X_1) +...+ coef_N*log(X_N)}, where \eqn{X_i ~ ChiSquare(df_i,\delta_i)}
#' are inedependent RVs, with \eqn{df_i} degrees of freedom and the noncentrality parameters \eqn{\delta_i >0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y = log(X)} with \eqn{X ~ ChiSquare(df,\delta)} is Poisson mixture
#' of CFs of the central log-transformed ChiSquare RVs of the form
#' \deqn{cf(t) = cf_LogRV_ChiSquareNC(t,df,\delta) = exp(-\delta/2) sum_{j=1}^Inf (\delta/2)^j/j! * cf_LogRV_ChiSquare(df+2*j),}
#' where \eqn{cf_LogRV_ChiSquare(df+j)} are the CFs of central log-transformed ChiSquare RVs with \eqn{df+j} degrees of freedom.
#' Hence, the characteristic function of \eqn{Y  = coef(1)*Y1 + ... + coef(N)*YN}
#' \deqn{cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t),}
#' where \eqn{cf_Yi(t)} is evaluated with the parameters \eqn{df_i} and \eqn{delta_i}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param df vector of the  degrees of freedom \code{df > 0}. If empty, default value is \code{df = 1}.
#' @param delta vector of the non-centrality parameters \code{delta > 0}. If empty, default value is \code{delta = 0}.
#' @param coef vector of the coefficients of the linear combination of the Beta distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#' @param tol tolerance factor for selecting the Poisson weights, i.e. such that \eqn{PoissProb > tol}.
#' If empty, default value is \code{tol = 1e-12}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent
#' LOG-TRANSFORMED non-central ChiSquare random variables.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution}.
#'
#' @family Continuous Probability Distribution
#' @family Non-central Probability Distribution
#'
#' @example R/Examples/example_cf_LogRV_ChiSquareNC.R
#'
#' @export
#'
cf_LogRV_ChiSquareNC <- function(t, df, delta, coef, niid, tol) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(df)) {
    df <- vector()
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

  if(length(df) == 0) {
    df <- 1
  }
  if(length(delta) == 0) {
    delta <- 0
  }
  if(length(coef) == 0) {
    coef <- 1
  }
  if(length(tol) == 0) {
    tol <- 1e-12
  }

  ## SET THE COMMON SIZE of the parameters
  l_max <- max(c(length(df), length(delta), length(coef)))
  if (l_max > 1) {
    if (length(df) == 1) {
      df <- rep(df, l_max)
    }
    if (length(delta) == 1) {
      delta <- rep(delta, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(
      coef, df, delta
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function of a linear combination of independent nc RVs
  szt <- dim(t)
  t <- c(t)
  cf <- rep(1, length(t))
  for(i in 1:length(coef)) {
    cf <- cf * cf_ncLogRVChiSquare(coef[i] * t, df[i], delta[i], tol)
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
cf_ncLogRVChiSquare <- function(t, df, delta, tol) {
  # cf_ncLogRVChiSquare Characteristic function of the distribution of the
  # non-central log-transformed ChiSquare RV with df degrees of freedom
  # and the non-centrality parameter delta > 0.

  f <- 0
  delta <- delta / 2
  if(delta == 0) {
    # Deal with the central distribution
    f <- cf_LogRV_ChiSquare(t, df)
  } else if(delta > 0) {
    # Sum the Poisson series of CFs of central iid ChiSquare RVs,
    # poisspdf(j,delta) .* cf_LogRV_ChiSquare(t,alpha+j,beta)
    j0 <- floor(delta / 2)
    p0 <- exp(-delta + j0 * log(delta) - log(gamma(j0 + 1)))
    f <- f + p0 * cf_LogRV_ChiSquare(t, df + 2 * j0)
    p <- p0
    j <- j0 - 1
    while(j >= 0 && p > tol) {
      p <- p * (j + 1) / delta
      f <- f + p * cf_LogRV_ChiSquare(t, df + 2 * j)
      j <- j - 1
    }
    p <- p0
    j <- j0+1
    i <- 0
    while(p > tol && i <= 5000) {
      p <- p * delta / j
      f <- f + p * cf_LogRV_ChiSquare(t, df + 2 * j)
      j <- j + 1
      i <- i + 1
    }
    if(i == 5000) {
      warning("No convergence.")
    }
    return(f)
  } else {
    stop('delta should be nonnegative.')
  }
}
