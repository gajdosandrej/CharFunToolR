#' @title Characteristic function of a linear combination of independent random variables
#' with the central FISHER-SNEDECOR F-distribution
#'
#' @description
#' \code{cf_FisherSnedecor(t, df1, df2, coef, niid, tol)} evaluates characteristic function
#' of the distribution of a linear combination of independent random variables
#' with the central FISHER-SNEDECOR F-distribution.
#'
#' That is, \code{cf_FisherSnedecor} evaluates the characteristic function \eqn{cf(t)}
#' of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ F(df1_i,df2_i)} are inedependent RVs,
#' with \eqn{df1_i} and \eqn{df2_i} degrees of freedom, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ F(df1,df2)} is
#' \deqn{cf(t) = U(df1/2, 1-df2/2, -1i*(df2/df1)*t),}
#' where \eqn{U(a,b,z)} denotes the confluent hypergeometric function of the second kind.
#'
#' Hence, the characteristic function of \eqn{Y  = coef_1*X_1 +...+ coef_N*X_N}
#' is \deqn{cf_Y(t) =  cf_1(coef_1*t) *...* cf_N(coef_N*t),} where \eqn{cf_i(t)}
#' is the characteristic function of \eqn{X_i ~ F(df1_i,df2_i)}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param df1 vector of the  degrees of freedom \code{df1 > 0}. If empty, default value is \code{df1 = 1}.
#' @param df2 vector of the  degrees of freedom \code{df2 > 0}. If empty, default value is \code{df2 = 1}.
#' @param coef vector of the coefficients of the linear combination of the log-transformed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#' @param tol tolerance factor for selecting the Poisson weights, i.e. such that \eqn{PoissProb > tol}.
#' If empty, default value is \code{tol = 1e-12}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent random variables
#' with the central FISHER-SNEDECOR F-distribution.
#'
#' @references
#' [1] PHILLIPS, P.C.B. The true characteristic function of the F distribution. Biometrika (1982), 261-264.
#'
#' [2] WITKOVSKY, V.: On the exact computation of the density and of the quantiles of linear combinations
#' of t and F random variables. Journal of Statistical Planning and Inference 94 (2001), 1–13.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/F-distribution}.
#'
#' @family Continuous Probability Distribution
#'
#' @example R/Examples/example_cf_FisherSnedecor.R
#'
#' @export
#'
cf_FisherSnedecor <- function(t, df1, df2, coef, niid, tol) {
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
  if(missing(tol)) {
    tol <- numeric()
  }

  if(length(df2) == 0 && length(df1) > 0) {
    df2 <- 1
  } else if (length(df2) == 0 && length(coef) > 0) {
    df2 <- 1
  } else if(any(df2 == 0) || length(df2) == 0) {
    df2 <- 1
  }

  if(length(df1) == 0 && length(coef) > 0) {
    df1 <- 1
  } else if (length(df1) == 0 && length(df2) > 0) {
    df1 <- 1
  }

  if(length(coef) == 0 && length(df2) > 0) {
    coef <- 1
  } else if (length(coef) == 0 && length(df1) > 0) {
    coef <- 1
  }

  if(length(tol) == 0) {
    tol <- 1e-6
  }
  reltol <- tol

  ## SET THE COMMON SIZE of the parameters
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

  ## Characteristic function of a linear combination of independent F RVs
  szt <- dim(t)
  t <- c(t)
  cf <- rep(1, length(t))
  for(i in 1:length(coef)) {
    cf_aux_re <- vector()
    cf_aux_im <- vector()
    for(j in 1:length(t)) {
      cf_aux_re <- c(cf_aux_re, integrate(function(x) {as.numeric(Re(funCF(df1[i], df2[i], coef[i] * t[j], (x / (1 - x)) ^ 2))) * (2 * x / (1 - x) ^ 3)}, 0, 1, rel.tol = reltol)$value)
      cf_aux_im <- c(cf_aux_im, integrate(function(x) {as.numeric(Im(funCF(df1[i], df2[i], coef[i] * t[j], (x / (1 - x)) ^ 2))) * (2 * x / (1 - x) ^ 3)}, 0, 1, rel.tol = reltol)$value)
      }
    cf <- cf * (cf_aux_re  + 1i*cf_aux_im)
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
funCF <- function(df1, df2, t, x) {
  # characteristic function CF of the F-distribution with df1 and df2
  # degrees of freedom at the real argument t.

  # The characteristic function of F-distribution with df1 and df2 degrees
  # of freedom, evaluated at the real t from (-inf,+inf), is defined by
  # CF(t) = (gamma(df1/2+df2/2)/gamma(df2/2)) * HypergeometricU(df1/2, 1-df2/2, -1i*(df2/df1)*t),
  # where HypergeometricU(a,b,z) denotes the confluent hypergeometric function of the second kind: U(a,b,z).

  # Here we use an integral representation of the hypergeometric function
  # U(a,b,z), defined for purly complex argument z as (VW2016):
  # U(a,b,z) = gamma(1-b)/gamma(a-b+1) * Integral_0^inf cfFun(a,b,x) dx.

  ## ALGORITHM FUNCF
  t <- c(t)
  nt <- length(t)
  o1 <- rep(1, nt)
  x <- t(Conj(x))
  nx <- length(x)
  o2 <- rep(1, nx)
  a <- df1 / 2
  b <- 1 - df2 / 2
  c <- (log(gamma(a - b + 1)) - log(gamma(1 - b)) - log(gamma(a)))
  z <- -(df2 / df1) * t
  f <- matrix(0, nt, nx)

  # z == 0
  id <- z == 0
  if(any(id > 0)) {
    f[id, ] <- exp(c + (a - 1) * log(x) + (b - a - 1) * log(1 + x))
  }

  # abs(z) >= 1
  id <- abs(z) >= 1
  if(any(id > 0)) {
    zi <- -1i / z[id]
    f[id,] <- exp(c + log(zi %*% t(o2)) + (a - 1) * log(zi %*% x) + (b - a - 1) * log(1 + zi %*% x) - o1[id] %*% x)
  }

  # z > 0 & z < 1
  id <- z>0 && z<1
  if(any(id > 0)) {
    f[id,] <- exp(c + log(-1i) + (a - 1) * log(-1i * o1[id] %*% x) + (b - a - 1) * log(1 - 1i * o1[id] %*% x) - (z[id]) %*% x)
  }

  # z < 0 & z > -1
  id <- z<0 && z>-1
  if(any(id > 0)) {
    f[id, ] <- exp(c + log(1i) + (a - 1) * log(1i * o1[id] %*% x) + (b - a - 1) * log(1 + 1i * o1[id] %*% x) + (z[id]) %*% x)
  }

  return(f)

}

## EXAMPLE
# df1 <- 5
# df2 <- 4
# t <- 1:5
# x <- seq(0, 1, length.out = 100)
# plotGraf(function(x)
#   funCF(df1, df2, t, x)[3,], x, title = '')
