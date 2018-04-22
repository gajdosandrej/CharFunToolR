#' @title Characteristic function of a linear combination of independent EXPONENTIAL random variables
#'
#' @description
#' \code{cf_Exponential(t, lambda)} evaluates characteristic function of a linear combination
#' (resp. convolution) of independent EXPONENTIAL random variables.
#'
#' That is, \code{cf_Exponential} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ EXP(\lambda_i)} are inedependent RVs,
#' with the rate parameters \eqn{\lambda_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y} is defined by
#' \deqn{cf(t) = Prod( \lambda_i / (\lambda_i - 1i*t) ).}
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Exponential_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param lambda vector of the 'rate' parameters \code{lambda > 0}. If empty, default value is \code{lambda = 1}.
#' @param coef vector of the coefficients of the linear combination
#' of the GAMMA random variables. If \code{coef} is scalar,
#' it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y +...+ Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent EXPONENTIAL random variables.
#'
#' @example R/Examples/example_cf_Exponential.R
#'
#' @export
cf_Exponential <- function(t, lambda, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(lambda)) {
    lambda <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  if (length(lambda) == 0) {
    lambda <- 1
  }
  if (length(coef) == 0) {
    coef <- 1
  }
  if (length(niid) == 0) {
    niid <- 1
  }

  # Equal size of the parameters
  if (length(coef) > 0 && length(lambda) == 1 && length(niid) == 0) {
    coef <- sort(coef)
    m <- length(coef)
    coef <- unique(coef)
    idx <- firstOccIdx(coef)
    lambda <- lambda * diff(idx, differences = m + 1)
  }

  l_max <- max(c(length(lambda), length(coef)))
  if (l_max > 1) {
    if (length(lambda) == 1) {
      lambda <- rep(lambda, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (any(lengths(list(coef, lambda)) < l_max)) {
      stop("Input size mismatch.")
    }
  }

  # Special treatment for linear combinations with large number of RVs
  szcoefs <- dim(as.matrix(coef))
  szcoefs <- szcoefs[1] * szcoefs[2]
  szt <- dim(as.matrix(t))
  sz <- szt[1] * szt[2]
  szcLimit <- ceiling(1e3 / (sz / (2 ^ 16)))
  idc <- 1:(as.integer(szcoefs / szcLimit) + 1)

  ## Characteristic function
  szt <- dim(t)
  t <- c(t)
  idx0 <- 1
  cf <- 1
  for (j in 1:idc[length(idc)]) {
    idx1 <- min(idc[j] * szcLimit, szcoefs)
    idx <- idx0:idx1
    idx0 <- idx1 + 1
    aux <- t %*% t(coef[idx] / lambda[idx])
    cf <- cf * apply(1/(1 - 1i * aux), 1, prod)
  }
  dim(cf) <- szt
  cf[t == 0] <- 1

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value")
    }
  }

  return(cf)

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

