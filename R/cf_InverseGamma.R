#' @title Characteristic function of a linear combination of independent INVERSE-GAMMA random variables
#'
#' @description
#' \code{cf_InverseGamma(t, alpha, beta, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of a linear combination (resp. convolution) of independent INVERSE-GAMMA random variables.
#'
#' That is, \code{cf_InvGamma} evaluates the characteristic function \eqn{cf(t)}
#' of  \eqn{Y =  sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ InvGamma(\alpha_i,\beta_i)}
#' are inedependent RVs, with the shape parameters \eqn{\alpha_i > 0} and
#' the  rate parameters \eqn{\beta_i > 0}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of Y is defined by
#' \deqn{cf(t) = Prod( 2 / gamma(\alpha(i)) * (-1i*\beta(i)*t).^(\alpha(i)/2) * besselk(\alpha(i),sqrt(-4i*\beta(i)*t)) ).}
#'
#' @family Continuous Probability Distribution
#'
#' @references
#' WITKOVSKY, V.: Computing the distribution of a linear combination
#' of inverted gamma variables, Kybernetika 37 (2001), 79-90.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha the shape parameter \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta the rate (1/scale) parameter \code{beta > 0}. If empty, default value is \code{beta = 1}.
#' @param coef  - vector of the coefficients of the linear combination
#' of the IGamma random variables. If coef is scalar,
#' it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y +...+ Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef_i * X_i}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @details
#' PARAMETRIZATION:\cr
#' As for the GAMMA distribution, also for the Inverse-Gamma distribution
#' there are three different possible parametrizations:\cr
#' i)   With a shape parameter \eqn{k} and a scale parameter \eqn{\theta}.\cr
#' ii)  With a shape parameter \eqn{\alpha = k} and an inverse scale parameter
#' \eqn{\beta = 1/\theta}, called a rate parameter.\cr
#' iii) With a shape parameter \eqn{k} and a mean parameter \eqn{\mu = k/\beta}.
#' In each of these three forms, both parameters are positive real numbers.
#'
#' Here, \code{cf_InverseGamma} implements the shape-rate parametrization with
#' parameters \eqn{\alpha} and \eqn{\beta}, respectively.
#'
#' If \eqn{X ~ IGamma(df/2,1/2)}(shape-rate parametrization), then \eqn{X} is identical
#' to \eqn{ICHI2(df)}, the Inverse-Chi-squared distribution with \eqn{df} degrees of freedom.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent INVERSE-GAMMA random variables.
#'
#' @note Ver.: 16-Sep-2018 18:28:11 (consistent with Matlab CharFunTool v1.3.0, 10-May-2017 18:11:50).
#'
#' @example R/Examples/example_cf_InverseGamma.R
#'
#' @export
#'
cf_InverseGamma <- function(t, alpha, beta, coef, niid) {
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
    niid <- numeric()
  }

  if (length(beta) == 0 && length(alpha) > 0) {
    beta <- 2
  } else if (length(beta) == 0 && length(coef) > 0) {
    beta <- 2
  } else if (any(beta == 0) || length(beta) == 0) {
    beta <- 2
  }

  if (length(alpha) == 0 && length(coef) > 0) {
    alpha <- 2
  } else if (length(alpha) == 0 && length(beta) > 0) {
    alpha <- 2
  }

  if (length(coef) == 0 && length(beta) > 0) {
    coef <- 1
  } else if (length(coef) == 0 && length(alpha) > 0) {
    coef <- 1
  }

  ## Equal size of the parameters
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
    aux <- ((-1i * t) %*% t((coef[idx] / beta[idx]))) ^ (1 / 2)
    aux0 <- matrix(0,dim(aux)[1],dim(aux)[2])
    aux1 <- 2*t(t(aux) ^ alpha[idx])
    aux2 <- rep(1,length(t)) %*% t(alpha[idx])
    aux3 <- 2 * aux
    for(k in 1:dim(aux)[1]) {
      for(l in 1:dim(aux)[2]) {
        aux0[k,l] <- tryCatch(Bessel::BesselK(aux3[k,l], aux2[k,l]), error = function(e) 0)
      }
    }
    aux <- aux1 * aux0
    aux <- t(t(aux) * (1 / gamma(alpha[idx])))
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
