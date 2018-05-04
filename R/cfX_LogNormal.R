#' @title Characteristic function of Lognormal distribution
#'
#' @description
#' \code{cfX_LogNormal(t, mu, sigma, tol)} computes the characteristic function \eqn{cf(t)}
#' of the Lognormal distribution with parameters \code{mu} (real) and \code{sigma > 0}, computed for real (vector) argument t, i.e.
#' \deqn{cf(t) = cfX_LogNormal(t,mu,sigma).}
#'
#' @importFrom stats integrate
#'
#' @details
#' In probability theory, a log-normal (or lognormal) distribution is a
#' continuous probability distribution of a random variable whose logarithm
#' is normally distributed. The lognormal distribution is defined for \eqn{x} in
#' \eqn{(0,+inf)} by its PDF/CDF/CF, as follows
#'  \eqn{pdf(x) = 1/(x*\sigma*sqrt(2*pi))*exp(-(ln(x)-\mu)^2/(2*\sigma^2))
#'  cdf(x) = 1/2+1/2*erf((ln(x)-\mu)/(sqrt(2)*\sigma))
#'  cf(t)  = sum_0^infinity{(it)^n/n!*exp(n*\mu + (n*\sigma)^2/2)}.}
#' As noted, this representation is asymptotically divergent but sufficient
#' for numerical purposes.
#'
#' \code{cfX_LogNormal} is based on the standard integral representation of the
#' characteristic function of the lognormal distribution, i.e.
#'  \eqn{cf(t) = Integral_0^inf exp(i*t*x) * PDF(x) dx}.
#' By using the half-space Fourier integral transformation we get
#'  \eqn{cf(t) = Integral_0^inf (i/t) * exp(-x) * PDF(i*x/t) dx}.
#' If we define the integrand as \eqn{funCF(t,x) = (i/t) * exp(-x) * PDF(i*x/t)},
#' then by using a stabilizing transformation from \eqn{[0,inf]} to \eqn{[0,1]}, we can
#' evaluate the CF by the following (well behaved) integral:
#'  \eqn{cf(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx}.
#'
#' \code{cfX_LogNormal} evaluates this integral by using the R built in
#' function \code{integrate()}, with precission specified by tolerance \code{tol} (default
#' value is \eqn{tol = 1e-6}).
#'
#' @family Continuous Probability Distribution
#'
#' @references
#' [1] WITKOVSKY, V.: On the exact computation of the density and
#' of the quantiles of linear combinations of t and F random variables.
#' Journal of Statistical Planning and Inference 94 (2001), 1-13.
#'
#' [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution of output quantity
#' in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' [3] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss distribution
#' based on numerical inversion of the compound empirical characteristic function
#' of frequency and severity. Working Paper. Insurance: Mathematics and Economics.
#'
#' [4] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for computing distributions
#' of collective risk models.  Working Paper. Journal of Statistical Software.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Log-normal_distribution}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param mu real, default value \code{mu = 0}.
#' @param sigma > 0, default value \code{sigma = 1}.
#' @param tol tolerance, default value \code{tol = 1e-6}.
#'
#' @return Characteristic function \eqn{cf(t)} of the Lognormal distribution.
#'
#' @example R/Examples/example_cfX_LogNormal.R
#'
#' @export
#'
cfX_LogNormal <- function(t,
                          mu = 0,
                          sigma = 1,
                          tol = 1e-6) {
  reltol <- tol

  szt <- dim(t)
  t <- c(t)

  cf <- seq(1, 1, length.out = length(t))
  cfIm <- seq(1, 1, length.out = length(t))
  cfRe <- seq(1, 1, length.out = length(t))
  id <- t != 0
  cfRe[id] <-
    unlist(lapply(
      t[id],
      FUN = function(t)
        integrate(function(x)
          Re(funCF(
            mu, sigma, t, (x / (1 - x)) ^ 2
          ) * (2 * x / (
            1 - x
          ) ^ 3)), 0, 1, rel.tol = reltol)$value
    ))
  cfIm[id] <-
    unlist(lapply(
      t[id],
      FUN = function(t)
        integrate(function(x)
          Im(funCF(
            mu, sigma, t, (x / (1 - x)) ^ 2
          ) * (2 * x / (
            1 - x
          ) ^ 3)), 0, 1, rel.tol = reltol)$value
    ))
  cf <- cfRe + 1i * cfIm

  dim(cf) <- szt

  return(cf)
}



# funCF Integrand function of the integral representation of the characteristic function CF
# of the Lognormal distribution with (scalar)  parameters mu (real) and sigma > 0 and the real (vector) argument t.

## EXAMPLE
# Plot the integrand functions for computing the CF(t)
# mu <- 0
# sigma <- 1
# t <- 1:5
# x <- seq(from = 0,to = 1,length.out = 100)
# f <- funCF(mu,sigma,t,x)
# plotGraf(function(t) funCF(mu,sigma,t,x), x, title = "")

# REFERENCES:
# WITKOVSKY V. (2016). On computing the characteristic functions
# of lognormal distribution and its applications. Working Paper.

funCF <- function(mu, sigma, t, x) {
  szt <- dim(t)
  t <- c(t)
  t  <- exp(mu) * t
  szx <- dim(x)
  x  <- c(x)

  ot <- seq(1, 1, length.out = length(t))
  dim(ot) <- szt
  ox <- seq(1, 1, length.out = length(x))
  dim(ox) <- szx

  funPDF <- function(x, s)
    exp(-0.5 * (log(x) / s) ^ 2) / (sqrt(2 * pi) * s * x)

  if (sigma >= 1 / 3) {
    t <- 1 / t
    f <- (1i * t * ox) * exp(-ot * x) * funPDF(1i * t * x, sigma)
  } else {
    # Set optimum limits for small and large abs(t)
    small <- 7 * sqrt(1 / sigma)
    large <- 25 * sqrt(1 / sigma)
    f <- ot * ox

    id <- abs(t) <= small
    if (any(id)) {
      f[id] <- exp(1i * t[id] * x) * funPDF(ot[id] * x * sigma)
    }

    id <- t > small & t <= large
    if (any(id)) {
      f[id] = exp(1i * ot[id] * x) * exp(-0.5 * (log(1 / t[id]) * x) / sigma) ^
        2 / (sqrt(2 * pi) * sigma * ot[id] * x)
    }

    id <- t < -small & t >= -large
    if (any(id)) {
      f[id] = exp(-1i * ot[id] * x) * exp(-0.5 * (log(-(1 / t[id]) * x) / sigma) ^
                                            2) / (sqrt(2 * pi) * sigma * ot[id] * x)
    }

    id <- abs(t) > large
    if (any(id)) {
      f[id] <- 0
    }

    return(f)

  }

}
