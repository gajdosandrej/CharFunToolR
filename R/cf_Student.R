#' @title Characteristic function of a linear combination
#' of independent STUDENT's t random variables
#'
#' @description
#' \code{cf_Student(t, df, mu, sigma, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
#' of a linear combination (resp. convolution) of independent (location and scale shifted) STUDENT's t random variables.
#'
#' That is, \code{cf_Student} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = sum_{i=1}^N coef_i * (mu_i + sigma_i * X_i)}, where \eqn{X_i ~ t(df_i)}
#' are inedependent (symmetric) t-distributed RVs, with \eqn{df_i > 0} degrees of freedom, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of the random variable \eqn{\mu + \sigma*X}, where \eqn{X ~ t(df)} is given by
#' \deqn{cf(t) = exp(1i*t*mu) * besselk(df/2,abs(sigma*t)*sqrt(df),1) * exp(-abs(sigma*t)*sqrt(df)) * (sqrt(df)*abs(aigma*t))^(df/2) / 2^(df/2-1)/gamma(df/2).}
#'
#' Hence, the characteristic function of \eqn{Y  = coef_1*(\mu_1+\sigma_1*X_1) +...+ coef_N*(\mu_N+\sigma_N*X_N)}
#' is \eqn{cf_Y(t) = exp(1i*\mu*t) * (cf_1(coef_1*\sigma_1*t) *...* cf_N(coef_N*\sigma_N*t))}, where \eqn{cf_i(t)}
#' is the characteristic function of \eqn{X_i ~ t(df_i)}.
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @references
#' WITKOVSKY V. (2016). Numerical inversion of a characteristic
#' function: An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Student's_t-distribution}.
#'
#' @param t  vector or array of real values, where the CF is evaluated.
#' @param df the degrees of freedom, \code{df > 0}. If empty, the default value is \code{df = 1}.
#' @param mu vector of location parameters, \code{mu} in Real. If empty, default value is \code{mu = 0}.
#' @param sigma vector of scale parameters, \code{sigma_i > 0}. If empty, default value is \code{sigma = 1}.
#' @param coef vector of the coefficients of the linear combination of the STUDENT's t random variables.
#' If \code{coef} is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient, such that \eqn{Z = Y + ... + Y} is
#' sum of \eqn{niid} random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef_i * X_i}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent STUDENT's t random variables.
#'
#' @example R/Examples/example_cf_Student.R
#'
#' @export
#'
cf_Student <- function(t, df, mu, sigma, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(df)) {
    df <- vector()
  }
  if (missing(mu)) {
    mu <- vector()
  }
  if (missing(sigma)) {
    sigma <- vector()
  }
  if (missing(coef)) {
    coef <- vector()
  }
  if (missing(niid)) {
    niid <- numeric()
  }

  if (length(df) == 0) {
    df <- 1
  }
  if (length(coef) == 0) {
    coef <- 1
  }
  if (length(mu) == 0) {
    mu <- 0
  }
  if (length(sigma) == 0) {
    sigma <- 1
  }
  if (length(niid) == 0) {
    niid <- 1
  }

  l_max <-
    max(c(length(coef), length(df), length(mu), length(sigma)))
  if (l_max > 1) {
    if (length(mu) == 1) {
      mu <- rep(mu, l_max)
    }
    if (length(sigma) == 1) {
      sigma <- rep(sigma, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (length(df) == 1) {
      df <- rep(df, l_max)
    }
    if ((any(lengths(list(
      coef, df, mu, sigma
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  szcoefs <- dim(as.matrix(coef))
  szcoefs <- szcoefs[1] * szcoefs[2]
  szt <- dim(as.matrix(t))
  sz <- szt[1] * szt[2]
  szcLimit <- ceiling(1e3 / (sz / (2 ^ 16)))
  idc <- 1:(as.integer(szcoefs / szcLimit) + 1)

  ## Characteristic function
  df2 <- df / 2
  szt <- dim(t)
  t <- c(t)
  o <- rep(1, length(t))
  idx0 <- 1

  for (j in 1:idc[length(idc)]) {
    idx1 <- min(idc[j] * szcLimit, szcoefs)
    idx <- idx0:idx1
    idx0 <- idx1 + 1
    aux <- abs(t) %*% t(abs(sqrt(df[idx]) * coef[idx] * sigma[idx]))
    aux0 <- matrix(0, dim(aux)[1], dim(aux)[2])
    o_df2 <- o %*% t(df2[idx])
    for(k in 1:dim(aux)[1]) {
      for(l in 1:dim(aux)[2]) {
        aux0[k,l] <- tryCatch(log(Bessel::BesselK(aux[k,l], o_df2[k,l], TRUE)), error = function(e) 0)
      }
    }
    aux <- - aux + t(df2[idx] * t(log(aux))) + aux0
    aux <- t(t(aux) + ((-log(2) * (df2[idx] - 1)) - log(gamma(df2[idx]))))
    aux <- 1i * t %*% t(coef[idx] * mu[idx]) + aux
    if(length(coef) > 1) {
      aux <- apply(aux, 1, sum)
    }
    cf <- exp(aux)
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
