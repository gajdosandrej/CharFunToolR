#' @title Characteristic function of a linear combination of independent GAMMA random variables
#'
#' @description
#' \code{cf_Gamma(t, alpha, beta, coef, niid)} evaluates the characteristic function
#' of a linear combination (resp. convolution) of independent GAMMA random variables.
#'
#' That is, \code{cf_Gamma} evaluates the characteristic function \eqn{cf(t)}
#' of  \eqn{Y =sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ Gamma(\alpha_i,\beta_i)}
#' are inedependent RVs, with the shape parameters \eqn{\alpha_i > 0}
#' and the rate parameters \eqn{\beta_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{Y} is defined by
#' \deqn{cf(t) = Prod( (1 - i*t*coef(i)/\beta(i))^(-\alpha(i)) ).}
#'
#' @family Continuous Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Gamma_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha the shape parameter \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta the rate (\eqn{1/scale}) parameter \code{beta > 0}. If empty, default value is \code{beta = 1}.
#' @param coef  - vector of the coefficients of the linear combination
#' of the GAMMA random variables. If coef is scalar,
#' it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#'  is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#'   is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @details
#' PARAMETRIZATION:\cr
#' Notice that there are three different parametrizations for GAMMA distribution in common use:\cr
#' i)   With a shape parameter \eqn{k} and a scale parameter \eqn{theta}.\cr
#' ii)  With a shape parameter \eqn{\alpha = k} and an inverse scale parameter \eqn{\beta = 1/\theta},
#' called a rate parameter.\cr
#' iii) With a shape parameter \eqn{k} and a mean parameter \eqn{\mu = k/\beta}.
#' In each of these three forms, both parameters are positive real numbers.\cr
#'
#' Here, \code{cf_Gamma} implements the shape-rate parametrization with parameters \code{alpha}
#' and \code{beta}, respectively.
#'
#' SPECIAL CASES:\cr
#' 1) If \eqn{X ~ Gamma(1,\lambda)} (shape-rate parametrization), then \eqn{X} has an exponential distribution
#' with rate parameter lambda.
#' 2) If \eqn{X ~ Gamma(df/2,1/2)}(shape-rate parametrization), then \eqn{X ~ ChiSquared(df)},
#' the chi-squared distribution with df degrees offreedom. Conversely,
#' if \eqn{Q ~ ChiSquared(df)} and c is a positive constant, then \eqn{cQ ~ Gamma(df/2,1/2c)}.
#' 3) If \eqn{X ~ Gamma(\alpha,\theta)} and \eqn{Y ~ Gamma(\beta,\theta)} are independently distributed,
#' then \eqn{X/(X + Y)} has a beta distribution with parameters \eqn{\alpha} and \eqn{\beta}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent GAMMA random variables.
#'
#' @note Ver.: 16-Sep-2018 18:26:34 (consistent with Matlab CharFunTool v1.3.0, 10-May-2017 18:11:50).
#'
#' @example R/Examples/example_cf_Gamma.R
#'
#' @export
#'
cf_Gamma <- function(t, alpha, beta, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(alpha)) {
    alpha <- vector()
  }
  if (missing(beta)) {
    beta <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- vector()
  }

  if (length(beta) == 0 && length(alpha) > 0) {
    beta <- 1
  } else if (length(beta) == 0 && length(coef) > 0) {
    beta <- 1
  } else if (any(beta == 0) || length(beta) == 0) {
    beta <- 1
  }

  if (length(alpha) == 0 && length(coef) > 0) {
    alpha <- 1
  } else if (length(alpha) == 0 && length(beta) > 0) {
    alpha <- 1
  } else if (any(alpha == 0) || length(alpha) == 0) {
    alpha <- 1
  }

  if (length(coef) == 0 && length(beta) > 0) {
    coef <- 1
  } else if (length(coef) == 0 && length(alpha) > 0) {
    coef <- 1
  }

  if (length(niid) == 0) {
    niid <- 1
  }

  ## Equal size of the parameters
  if (length(coef) > 0 &&
      length(alpha) == 1 && length(beta) == 1 && length(niid) == 0) {
    coef <- sort(coef)
    m <- length(coef)
    coef <- unique(coef)
    idx <- firstOccIdx(coef)
    lambda <- lambda * diff(idx, differences = m + 1)
  }

  l_max <- max(c(length(alpha), length(beta), length(coef)))
  if (l_max > 1) {
    if (length(alpha) == 1) {
      alpha <- rep(alpha, l_max)
    }
    if (length(beta) == 1) {
      beta <- rep(beta, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (any(lengths(list(coef, alpha, beta)) < l_max)) {
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
    aux <- t %*% t(coef[idx] / beta[idx])
    aux <- t(t((1 - 1i * aux)) ^ (-alpha[idx]))
    cf <- cf * apply(aux, 1, prod)
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
