#' @title Characteristic function  of a linear combination of independent CHI-SQUARE random variables
#'
#' @description
#' \code{cf_ChiSquare(t, df, ncp, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#'  of a linear combination (resp.convolution)
#'  of independent (possibly non-central) CHI-SQUARE random variables.
#'
#' That is, \code{cf_ChiSquare} evaluates the characteristic function \eqn{cf(t)}
#' of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ ChiSquare(df_i,ncp_i)}
#' are inedependent RVs, with \eqn{df_i > 0} degrees of freedom the 'non-centrality'
#' parameters \eqn{ncp_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y} is defined by
#' \deqn{cf(t) = Prod ( (1-2*i*t*coef(i))^(-df(i)/2) * exp((i*t*ncp(i))/(1-2*i*t*coef(i))) ).}
#'
#' @family Continuous Probability Distribution
#'
#' @references
#' IMHOF J. (1961): Computing the distribution of quadratic forms in normal variables. Biometrika 48, 419-426.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Chi-squared_distribution},\cr
#' \url{https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param df the degrees of freedom parameter \code{df > 0}. If empty, default value is \code{df = 1}.
#' @param ncp the non-centrality parameter \code{ncp > 0}. If empty, default value is \code{ncp = 0}.
#' @param coef vector of the coefficients of the linear combination
#' of the chi-squared random variables. If coef is scalar,
#' it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent CHI-SQUARE random variables.
#'
#' @note Ver.: 16-Sep-2018 18:16:32 (consistent with Matlab CharFunTool v1.3.0, 10-May-2017 18:11:50).
#'
#' @example R/Examples/example_cf_ChiSquare.R
#'
#' @export
#'
cf_ChiSquare <- function(t, df, ncp, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(df)) {
    df <- vector()
  }
  if (missing(ncp)) {
    ncp <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- vector()
  }

  isNoncentral <- FALSE

  if (length(ncp) == 0 && length(df) > 0) {
    isNoncentral <- FALSE
    ncp <- 0
  } else if (length(ncp) == 0 && length(coef) > 0) {
    isNoncentral <- FALSE
    ncp <- 0
  } else if(any(ncp==0) || length(ncp) == 0) {
    isNoncentral <- FALSE
    ncp <- 0
    } else {
    isNoncentral <- TRUE
  }

  if (length(df) == 0 && length(coef) > 0) {
    df <- 1
  } else if (length(df) == 0 && length(ncp) > 0) {
    df <- 1
  }

  if (length(coef) == 0 && length(ncp) > 0) {
    coef <- 1
  } else if (length(coef) == 0 && length(df) > 0) {
    coef <- 1
  }

  if (length(niid) == 0) {
    niid <- 1
  }

  firstOccIdx <- function(x) {
    x_uniqe <- sort(unique(x))
    indices <- vector()
    for (i in 1:length(x_uniqe)) {
      idx <- 1
      for (j in 1:length(x)) {
        if (x[j] == x_uniqe[i]) {
          idx <- j
          indices <- c(indices, idx)
          break
        }
      }
    }
    return(indices)
  }

  ## Find the unique coefficients and their multiplicities
  if (length(coef) > 0 && length(df) == 1 && isNoncentral == FALSE) {
    coef <- sort(coef)
    coef_orig_sort <- coef
    m <- length(coef)
    coef <- unique(coef)
    idx <- firstOccIdx(coef_orig_sort)
    df <- df * diff(c(idx,m+1))
  }

  # Check/set equal dimensions for the vectors coef, df, and ncp
  l_max <- max(c(length(df), length(ncp), length(coef)))
  if (l_max > 1) {
    if (length(df) == 1) {
      df <- rep(df, l_max)
    }
    if (length(ncp) == 1) {
      ncp <- rep(ncp, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (any(lengths(list(coef, df, ncp)) < l_max)) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function
  szt <- dim(t)
  t <- c(t)
  cf <- 1
  aux <- t %*% t(coef)
  #aux <- c(t %*% t(coef))
  dim_aux <- dim(t %*% t(coef))
  if (isNoncentral) {
    aux <- (t(t(1 - 2i * aux) ^ (-df / 2))) * exp((aux * (1i * ncp)) / (1 - 2i *
                                                                      aux))
  } else {
    aux <- t(t(1 - 2i * aux) ^ (-df / 2))
  }
  if(length(dim_aux)>0) {
    cf <- cf * apply(aux, 1, prod)

  } else {
    cf <- cf * aux
  }
  dim(cf) <- szt

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value")
    }
  }

  return(cf)

}
