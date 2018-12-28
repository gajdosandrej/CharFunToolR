#' @title Characteristic function of a linear combination
#' of independent LOG-TRANSFORMED CHI-SQUARE random variables
#'
#' @description
#' \code{cf_LogRV_ChiSquare(t, df, coef, niid)} evaluates characteristic function of a linear combination
#' (resp. convolution) of independent LOG-TRANSFORMED CHI-SQUARE random variables (RVs)
#' \eqn{log(X)}, where \eqn{X ~ ChiSquare(df)} is central CHI-SQUARE distributed RV with df degrees of freedom.
#'
#' That is, \code{cf_LogRV_ChiSquare} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = coef_1*log(X_1) +...+ coef_N*log(X_N)}, where \eqn{X_i ~ ChiSquare(df_i)}, with \eqn{df_i > 0}
#' degrees of freedom, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y = log(X)}, with \eqn{X ~ ChiSquare(df)} is defined by
#' \deqn{cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)).}
#' That is, the characteristic function can be derived from expression for the r-th moment of \eqn{X},
#' \eqn{E(X^r)} by using \eqn{(1i*t)} instead of \eqn{r}.
#' In particular, the characteristic function of \eqn{Y = log(X)} is
#' \deqn{cf_Y(t) = 2^(1i*t) * gamma(df/2 + 1i*t) / gamma(df/2).}
#'
#'Hence,the characteristic function of \eqn{Y  = coef_1*X_1 +...+ coef_N*X_N}
#'is  \deqn{cf_Y(t) =  cf_1(coef_1*t) *...* cf_N(coef_N*t),}
#'where \eqn{cf_i(t)} is the characteristic function of the ChiSquare distribution with \eqn{df_i} degrees of freedom.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param df vector of the  degrees of freedom \code{df > 0}. If empty, default value is \code{df = 1}.
#' @param coef vector of the coefficients of the linear combination of the logGamma random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent LOG-TRANSFORMED CHI-SQUARE random variables.
#'
#' @references
#' [1] PHILLIPS, P.C.B. The true characteristic function of the F distribution. Biometrika (1982), 261-264.
#'
#' [2] WITKOVSKY, V.: On the exact computation of the density and of the quantiles of linear combinations
#' of t and F random variables. Journal of Statistical Planning and Inference 94 (2001), 1–13.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Chi-squared_distribution}.
#'
#' @family Continuous Probability Distribution
#'
#' @note Ver.: 16-Sep-2018 18:34:40 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cf_LogRV_ChiSquare.R
#'
#' @export
#'
cf_LogRV_ChiSquare <- function(t, df, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(df)) {
    df <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- numeric()
  }

  ## SET the default values
  # if(length(df) == 0 && length(coef) > 0) {
  #   df <- 1
  # }
  # if(length(coef) == 0 && length(df) > 0) {
  #   coef <- 1
  # }
  # if(length(niid) == 0) {
  #   niid <- 1
  # }

  if(length(df) == 0) {
    df <- 1
  }
  if(length(coef) == 0) {
    coef <- 1
  }
  if(length(niid) == 0) {
    niid <- 1
  }

  ## Check size of the parameters
  l_max <- max(c(length(df), length(coef)))
  if (l_max > 1) {
    if (length(df) == 1) {
      df <- rep(df, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(
      coef, df
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function of linear combination
  szt <- dim(t)
  t   <- c(t)
  aux <- 1i * t %*% t(coef)
  aux <- GammaLog(t(t(aux) + df / 2)) - rep(1, length(t)) %*% t(GammaLog(df / 2))
  aux <- aux + 1i * t * log(2)
  cf  <- apply(exp(aux), 1, prod)
  dim(cf)  <- szt
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
