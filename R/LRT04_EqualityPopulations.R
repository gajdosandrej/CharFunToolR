#' @title
#' p-value of the log-transformed LRT statistic and/or its null
#' distribution CF/PDF/CDF
#'
#' @description
#' \code{LRT04_EqualityPopulations(W, n, p, q, options)} computes \eqn{p}-value of the log-transformed LRT
#' statistic \eqn{W = -log(\Lambda)}, for testing the null hypothesis of equality
#' of \eqn{q} (\eqn{q > 1}) \eqn{p}-dimensional normal populations, and/or its null distribution CF/PDF/CDF.
#' This is based on BALANCED samples of size n for each population!
#'
#' This is an EXPERIMENTAL version. Correctness should be checked again!
#'
#' @param W observed value of the minus log-transformed LRT statistic
#' \eqn{W = -log(\Lambda)}. If empty, the  algorithm evaluates the
#' CF/PDF/CDF and the quantiles of the null distribution of \eqn{W}.
#' @param n sample size, \eqn{n > min(p+q-1)}.
#' @param p common dimension of the vectors \eqn{X_k, k = 1,...q}.
#' @param q number of normal populations, \eqn{q > 1}.
#' @param options option structure, for more details see \code{\link{cf2DistGP}}. Moreover, \cr
#' \code{x} set vector of values where PDF/CDF is evaluated, \cr
#' \code{prob} set vector of probabilities for the quantiles, \cr
#' \code{coef}  set arbitrary multiplicator of the argument \code{t}
#' of the characteristic function. If empty, default value is \eqn{-n/2}
#' (standard value for minus log-transform of LRT). Possible
#' alternative is e.g. \code{coef = -1}, leading to \eqn{W = -(2/n)*log(LRT)}.
#'
#' @details
#' In particular, let \eqn{X_k ~ N_p(\mu_k,\Sigma_k)}, for \eqn{k = 1,...,q}. We want to
#' test the hypothesis that the \eqn{q} normal populations are equally
#' distributed. That is, we want to test that the
#' mean vectors \eqn{\mu_k} are equal for all \eqn{k = 1,...,q}, as well as the
#' covariance matrices \eqn{\Sigma_k} are equal for all \eqn{k = 1,...,q}. Then,
#' the null hypothesis is given as \eqn{H0: \mu_1 = ... = \mu_q} & \eqn{\Sigma_1 = ... = \Sigma_k}.
#' Here, the null hypothesis H0 and the LRT statistic can be decomposed:
#' \eqn{\Lambda = \Lambda_Means * \Lambda_Covariances}
#' where (first) \eqn{\Lambda_Covariances} represents the LRT for testing equality of
#' covariance matrices of given \eqn{q} normal populations, and (second)
#' \eqn{\Lambda_Means} represents (conditionally) the LRT for testing equality of
#' means of given \eqn{q} normal populations.
#'
#' Under null hypothesis, distributions of \eqn{\Lambda_Covariances} and
#' \eqn{\Lambda_Means} are independent, and the distribution of the test statistic
#' \eqn{\Lambda} is \eqn{Lambda ~  \Lambda_Means * \Lambda_Covariances},
#' \eqn{~  (prod_{k=1}^q prod_{j=1}^{p} (B_{jk})^{n/2})* (prod_{j=1}^{p} (B_j)^{n*q/2})}
#' where the \eqn{B_{jk}} and \eqn{B_j} are mutually independent beta distributed
#' random variables. Here we assume that \eqn{n} is equal sample size for each
#' sample, \eqn{k = 1,...,q}, \eqn{n > p}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic Lambda, say \eqn{W = -log(\Lambda)} is given by
#' \eqn{cf = function(t) {cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q))}}
#' . * cf_LogRV_Beta(-(n*q/2)*t, ((n-1)*q-i+1)/2, (q-1)/2),
#' where \eqn{i = (1, 2, ..., p)}, \eqn{k = (1*o,...,q*o)} with \eqn{p}-dimensional vector
#' of ones \eqn{o = (1,...,1)}  and \eqn{j = (j_1,...,j_q)} with \eqn{j_k = 1:p}.
#'
#' @return
#' \eqn{p}-value of the log-transformed LRT statistic, \eqn{W = -log(\Lambda)}
#' and/or its null distribution CF/PDF/CDF.
#'
#' @references
#' [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate Statistical Analysis.
#' New York: Wiley, 3rd Ed., 2003.
#'
#' [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
#' near-exact distribution theory for the most common likelihood ratio
#' test statistics used in Multivariate Analysis. Test, 2011, 20.1:180-203.
#'
#' [3] WITKOVSKY, Viktor. Exact distribution of selected multivariate test
#' criteria by numerical inversion of their characteristic functions. \cr
#' arXiv preprint arXiv:1801.02248, 2018.
#'
#' @family Likelihood Ratio Test
#'
#' @note Ver.: 16-Sep-2018 21:10:45 (consistent with Matlab CharFunTool v1.3.0, 20-Jan-2018 12:43:15).
#'
#' @example R/Examples/example_LRT04_EqualityPopulations.R
#'
#' @export
#'
LRT04_EqualityPopulations <- function(W, n, p, q, options) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(n) || missing(p) || missing(q)) {
    stop("Enter input parameters n, p, q.")
  }

  if (missing(options)) {
    options <- list()
  }

  if (is.null(options$x)) {
    options$x <- vector()
  }

  if (is.null(options$prob)) {
    options$prob <- vector()
  }

  if (is.null(options$coef)) {
    options$coef <- vector()
  }

  if (is.null(options$xMin)) {
    options$xMin <- 0
  }

  if (n <= p) {
    stop("Sample size n is too small.")
  }

  ind <- Conj(1:p)
  alpha <- rep(0, p * (q + 1))
  beta <- rep(0, p * (q + 1))
  coef <- rep(0, p * (q + 1))
  # Parameters of Beta distributions used for LRT_Means
  alpha[1:p] <- ((n - 1) * q - ind + 1) / 2
  beta[1:p] <- (q - 1) / 2
  coef[1:p] <- -n * q / 2        # CHECK THIS!!!
  # Parameters of Beta distributions used for LRT_Covariances
  coef[p + (1:(p * q))] <- -n / 2
  # CHECK THIS!!!
  ind <- p
  for (k in 1:q) {
    for (j in 1:p) {
      ind <- ind + 1
      alpha[ind] <- (n - j) / 2
      beta[ind] <- (j * (q - 1) + 2 * k - 1 - q) / (2 * q)
    }
  }
  cf <- function(t) {
    cf_LogRV_Beta(t, alpha, beta, coef)
  }

  # Evaluate the p-value and PDF/CDF/QF of the log-transformed LRT statistic.
  if (!missing(W) && length(W) > 0) {
    # % P-VALUE
    options$isPlot <- FALSE
    result <- cf2DistGP(cf = cf, x = W, options = options)
    pval <- 1 - result$cdf
  } else {
    # DISTRIBUTION of Lambda PDF/CDF
    pval <- vector()
    result <-
      cf2DistGP(
        cf = cf,
        x = options$x,
        prob = options$prob,
        options = options
      )
  }

  result$alpha <- alpha
  result$beta <- beta

  return(list(pval, result))

}
