#' @title Characteristic function of a linear combination of independent
#'  LOG-TRANSFORMED BETA random variables
#'
#' @description
#' \code{cf_LogRV_Beta(t, alpha, beta, coef, niid)} evaluates characteristic function of a linear combination (resp.
#' convolution) of independent LOG-TRANSFORMED BETA random variables (RVs)
#' \eqn{log(X)}, where \eqn{X ~ BETA(\alpha,\beta)}, and \eqn{\alpha > 0} and \eqn{\beta > 0}
#' represent the shape parameters of the BETA distribution.
#'
#' That is, \code{cf_LogRV_Beta} evaluates the characteristic function \eqn{cf(t)} of  \eqn{Y =}
#' \eqn{coef_i*log(X_1) +...+ coef_N*log(X_N)}, where \eqn{X_i ~ Beta(\alpha_i,\beta_i)}
#' are inedependent RVs, with the shape parameters \eqn{\alpha_i > 0} and \eqn{\beta_i > 0}, for \eqn{i = 1,...,N}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha vector of the 'shape' parameters \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta vector of the 'shape' parameters \code{beta > 0}. If empty, default value is \code{beta = 1}.
#' @param coef vector of the coefficients of the linear combination of the
#' log-transformed random variables. If \code{coef} is scalar, it is
#' assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y} is
#' sum of \code{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N}
#' \eqn{coef(i) * log(X_i)} is independently and identically distributed
#' random variable. If empty, default value is \code{niid = 1}.
#'
#' @details
#' The characteristic function of \eqn{Y = log(X)}, with \eqn{X ~ Beta(\alpha,\beta)} is
#' defined by \eqn{cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t))}.
#' That is, the characteristic function can be derived from expression for
#' the \eqn{r}-th moment of \eqn{X}, E(X^r) by using \eqn{(1i*t)} instead of \eqn{r}. In
#' particular, the characteristic function of \eqn{Y = log(X)}, with \eqn{X ~}
#' \eqn{Beta(\alpha,\beta)} is defined by
#' \eqn{cf_Y(t) = gamma(\alpha + 1i*t) / gamma(\alpha) .* ...}
#' \eqn{gamma(\alpha + \beta) / gamma(\alpha + \beta + 1i*t)}.
#' Hence,the characteristic function of \eqn{Y  = coef(1)*Y1 + ... + coef(N)*YN}
#' is \eqn{cf_Y(t) = cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t)}, where \eqn{cf_Yi(t)}
#' is evaluated with the parameters \code{alpha(i)} and \code{beta(i)}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent
#'  LOG-TRANSFORMED BETA random variables.
#'
#' @references
#' GLASER, R.E. (1976). The ratio of the geometric mean to the arithmetic
#' mean for a random sample from a gamma distribution. \emph{Journal of the American Statistical Association}, 71(354), 480-487.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}.
#'
#' @family Continuous Probability Distribution
#'
#' @note Ver.: 16-Sep-2018 18:29:18 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
#'
#' @example R/Examples/example_cf_LogRV_Beta.R
#'
#' @export
#'
cf_LogRV_Beta <- function(t, alpha, beta, coef, niid) {
  ##CHECK THE INPUT PARAMETERS
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

  ##
  if (length(beta) == 0 && length(alpha) > 0) {
    beta <- 1
  } else if (length(beta) == 0 && length(coef) > 0) {
    beta <- 1
  }
  if (length(alpha) == 0 && length(coef) > 0) {
    alpha <- 1
  } else if (length(alpha) == 0 && length(beta) > 0) {
    alpha <- 1
  }
  if (length(coef) == 0 && length(beta) > 0) {
    coef <- 1
  } else if (length(coef) == 0 && length(alpha) > 0) {
    coef <- 1
  }

  ## Check size of the parameters
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
    if ((any(lengths(list(
      coef, alpha, beta
    )) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function of a linear combination
  szt <- dim(t)
  t <- c(t)
  aux <- 1i * (t %*% t(coef))
  if (length(t) == 1) {
    aux <- GammaLog((t(aux + alpha))) - GammaLog((t(aux + (alpha + beta))))
    aux <- aux + (GammaLog(alpha + beta) - GammaLog(alpha))
    cf <- apply(t(exp(aux)), 1, prod)
  } else {
    aux <-
      GammaLog(sweep(aux, 2, alpha, "+")) - GammaLog(sweep(aux, 2, (alpha + beta), "+"))
    aux <-
      aux + (rep(1, length(t)) %*% t((
        GammaLog(alpha + beta) - GammaLog(alpha)
      )))
    cf <- apply(exp(aux), 1, prod)
  }

  dim(cf) <- szt
  cf[t == 0] <- 1

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value.")
    }
  }
  return(cf)
}
