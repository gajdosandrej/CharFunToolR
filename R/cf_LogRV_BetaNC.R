#' @title Characteristic function of a linear combination
#' of independent LOG-TRANSFORMED non-central BETA random variables
#'
#' @description
#' \code{cf_LogRV_BetaNC(t, alpha, beta, delta, coef, niid, tol, type)} evaluates characteristic function
#' of a linear combination (resp. convolution) of independent LOG-TRANSFORMED non-central BETA random variables,
#' (Type I and Type II), with their distributions specified
#' by the parameters \eqn{\alpha_i, \beta_i}, and the noncentrality parameters \eqn{\delta_i}.
#'
#' The noncentral beta distribution has two types.
#' The Type I is the distribution  of the random variable \eqn{B1 = X1/(X1+X2)}, \eqn{X1 ~ Gamma(\alpha, \gamma, \delta)}
#' and \eqn{X2 ~ Gamma(\beta, \gamma)}. The Type II noncentral beta distribution is the distribution
#' of the ratio random variable \eqn{B2 = X1/(X1+X2)}, where \eqn{X1 ~ Gamma(\alpha, \gamma)} and
#'\eqn{X2 ~ Gamma(\beta, \gamma, \delta)}.
#'
#' That is, \code{cf_LogRV_BetaNC} evaluates the characteristic function \eqn{cf(t)}
#' of \eqn{Y = coef_i*log(X_1) +...+ coef_N*log(X_N)}, where \eqn{X_i ~ Beta(\alpha_i,\beta_i,\delta_i)}
#' are inedependent RVs, with the shape parameters \eqn{\alpha_i > 0}, \eqn{\beta_i > 0},
#' and the noncentrality parameters \eqn{\delta_i > 0}, for \eqn{i = 1,...,N}.
#'
#' For Type I noncentral beta distribution, the characteristic
#' function of \eqn{Y = log(X)} with \eqn{X ~ Beta(\alpha, \beta, \delta)} is Poisson mixture
#' of CFs of the central log-transformed Beta RVs of the form
#' \deqn{cf(t) = cf_LogRV_BetaNC(t, \alpha, \beta, \delta) = exp(-\delta/2) * sum_{j=1}^Inf (\delta/2)^j/j! * cf_LogRV_Beta(\alpha+j, \beta)}
#' where \eqn{cf_LogRV_Beta(\alpha+j, \beta)} are the CFs of central log-transformed Beta RVs
#' with parameters \eqn{\alpha+j} and \eqn{\beta}. Alternatively,
#' \deqn{cf(t) = cf_LogRV_BetaNC(t, \alpha, \beta, \delta) = Gamma(\alpha+1i*t)/Gamma(\alpha) * Gamma(\alpha+\beta)/Gamma(\alpha+\beta+1i*t) * exp(-\delta/2) * 2F2(\alpha+\beta, \alpha+1i*t; \alpha, \alpha+\beta+1i*t; \delta/2)},
#' where \eqn{2F2(a,b,;c,d;z)} is hypergeometric function.
#'
#' For Type II noncentral beta distribution, the characteristic function
#' of \eqn{Y = log(X)} with \eqn{X ~ Beta(\alpha, \beta, \delta)} is Poisson mixture of CFs
#' of the central log-transformed Beta RVs of the form
#' \deqn{cf(t) = cf_LogRV_BetaNC(t, \alpha, \beta, \delta) = exp(-\delta/2) * sum_{j=1}^Inf (\delta/2)^j/j! * cf_LogRV_Beta(\alpha, \beta+j)}
#' where \eqn{cf_LogRV_Beta(\alpha, \beta+j)} are the CFs of central log-transformed Beta RVs with parameters \eqn{\alpha} and \eqn{\beta+j}.
#' Alternatively,
#' \deqn{cf(t) = cf_LogRV_BetaNC(t, \alpha, \beta, \delta) = Gamma(\alpha+1i*t)/Gamma(\alpha) * Gamma(\alpha+\beta)/Gamma(\alpha+ \beta+1i*t) * 1F1(1i*t;1i*t+(\alpha+\beta); -\delta/2) = cf_LogRV_Beta(t, \alpha, \beta) * 1F1(1i*t;1i*t+(\alpha+\beta);-\delta/2)}
#' where \eqn{1F1(a;b;z)} is hypergeometric function.
#'
#' Hence,the characteristic function of \eqn{Y  = coef(1)*Y1 + ... + coef(N)*YN}
#' is  \eqn{cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t)}, where \eqn{cf_Yi(t)} is evaluated
#' with the parameters \eqn{\alpha(i)} and \eqn{\beta(i)}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha vector of the 'shape' parameters \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta vector of the 'shape' parameters \code{beta > 0}. If empty, default value is \code{beta = 1}.
#' @param delta vector of the non-centrality parameters \code{delta > 0}. If empty, default value is \code{delta = 0}.
#' @param coef vector of the coefficients of the linear combination of the Beta distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#' @param tol tolerance factor for selecting the Poisson weights, i.e. such that \eqn{PoissProb > tol}.
#' If empty, default value is \code{tol = 1e-12}.
#' @param type indicator of the type of the noncentral distribution
#' (Type I = 1 or Type II = 2). If empty, default value is \code{type = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination
#' of independent LOG-TRANSFORMED non-central BETA random variables.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Noncentral_beta_distribution}.
#'
#' @family Continuous Probability Distribution
#' @family Non-central Probability Distribution
#'
#' @note Ver.: 20-Sep-2018 00:17:22 (consistent with Matlab CharFunTool v1.3.0, 06-Jul-2018 15:10:17).
#'
#' @example R/Examples/example_cf_LogRV_BetaNC.R
#'
#' @export
#'
cf_LogRV_BetaNC <- function(t, alpha, beta, delta, coef, niid, tol, type) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(alpha)) {
    alpha <- vector()
  }
  if(missing(beta)) {
    beta <- vector()
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
  if(missing(type)) {
          type <- numeric()
  }

  if(length(alpha) == 0) {
    alpha <- 1
  }
  if(length(beta) == 0) {
    beta <- 1
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
  if(length(type) == 0) {
          type <- 1
  }

  ## SET THE COMMON SIZE of the parameters
  l_max <- max(c(length(alpha), length(beta), length(delta), length(coef)))
  if (l_max > 1) {
    if (length(alpha) == 1) {
      alpha <- rep(alpha, l_max)
    }
    if (length(beta) == 1) {
      beta <- rep(beta, l_max)
    }
    if (length(delta) == 1) {
      delta <- rep(delta, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if ((any(lengths(list(
      coef, alpha, beta, delta
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function of a linear combination of independent nc F RVs
  szt <- dim(t)
  t <- c(t)
  cf <- rep(1, length(t))
  for(i in 1:length(coef)) {
    cf <- cf * cf_ncLogRVBeta(coef[i] * t, alpha[i], beta[i], delta[i], tol, type)
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
cf_ncLogRVBeta <- function(t, alpha, beta, delta, tol, type) {
  # cf_ncLogRVBeta Characteristic function of the distribution of the
  # non-central log-transformed Beta RV with the parameters alpha and beta,
  # and the non-centrality parameter delta > 0.

  f <- 0
  delta <- delta / 2
  if(delta == 0) {
    # Deal with the central distribution
    f <- cf_LogRV_Beta(t, alpha, beta)
  } else if(delta > 0) {
          if(type == 1) { # Type I noncentral beta distribution
                  # Sum the Poisson series of CFs of central iid Beta RVs,
                  # poisspdf(j,delta) .* cf_LogRV_Beta(t,alpha+j,beta)
                  j0 <- floor(delta / 2)
                  p0 <- exp(-delta + j0 * log(delta) - log(gamma(j0 + 1)))
                  f <- f + p0 * cf_LogRV_Beta(t, alpha + j0, beta)
                  p <- p0
                  j <- j0 - 1
                  while(j >= 0 && p > tol) {
                          p <- p * (j + 1) / delta
                          f <- f + p * cf_LogRV_Beta(t, alpha + j, beta)
                          j <- j - 1
                  }
                  p <- p0
                  j <- j0 + 1
                  i <- 0
                  while(p > tol && i <= 5000) {
                          p <- p * delta / j
                          f <- f + p * cf_LogRV_Beta(t, alpha + j, beta)
                          j <- j + 1
                          i <- i + 1
                  }
                  if(i == 5000) {
                          warning("No convergence.")
                  }
                  return(f)

          } else { # Type I noncentral beta distribution
                  # Sum the Poisson series of CFs of central iid Beta RVs,
                  # poisspdf(j,delta) .* cf_LogRV_Beta(t,alpha,beta+j)
                  j0 <- floor(delta / 2)
                  p0 <- exp(-delta + j0 * log(delta) - log(gamma(j0 + 1)))
                  f <- f + p0 * cf_LogRV_Beta(t, alpha, beta + j0)
                  p <- p0
                  j <- j0 - 1
                  while(j >= 0 && p > tol) {
                          p <- p * (j + 1) / delta
                          f <- f + p * cf_LogRV_Beta(t, alpha, beta + j)
                          j <- j - 1
                  }
                  p <- p0
                  j <- j0 + 1
                  i <- 0
                  while(p > tol && i <= 5000) {
                          p <- p * delta / j
                          f <- f + p * cf_LogRV_Beta(t, alpha, beta + j)
                          j <- j + 1
                          i <- i + 1
                  }
                  if(i == 5000) {
                          warning("No convergence.")
                  }
                  return(f)

          }

  } else {
    stop("delta should be nonnegative.")
  }

}
