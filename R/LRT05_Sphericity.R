#' @title
#' p-value of the log-transformed LRT statistic and/or its null
#' distribution CF/PDF/CDF
#'
#' @description
#' \code{LRT05_Sphericity(W, n, p, options)} computes \eqn{p}-value of the log-transformed LRT
#' statistic \eqn{W = -log(\Lambda)}, for testing the null hypothesis of equality
#' of \eqn{q} (\eqn{q > 1}) \eqn{p}-dimensional normal populations, and/or its null distribution CF/PDF/CDF.
#'
#' @param W observed value of the minus log-transformed LRT statistic
#' \eqn{W = -log(\Lambda)}. If empty, the  algorithm evaluates the
#' CF/PDF/CDF and the quantiles of the null distribution of \eqn{W}.
#' @param n sample size, \eqn{n > min(p+q-1)}.
#' @param p common dimension of the vectors \eqn{X_k, k = 1,...q}.
#' @param options option structure, for more details see \code{\link{cf2DistGP}}. Moreover, \cr
#' \code{x} set vector of values where PDF/CDF is evaluated, \cr
#' \code{prob} set vector of probabilities for the quantiles, \cr
#' \code{coef}  set arbitrary multiplicator of the argument \code{t}
#' of the characteristic function. If empty, default value is \eqn{-n/2}
#' (standard value for minus log-transform of LRT). Possible
#' alternative is e.g. \code{coef = -1}, leading to \eqn{W = -(2/n)*log(LRT)}.
#'
#' @details
#' In particular, let \eqn{X ~ N_p(\mu,\Sigma)}, then, the null hypothesis is given as
#' \eqn{H0: \Sigma = \sigma^2 I_p} (\eqn{\sigma^2} unspecified).
#' Here, the LRT test statistic is given by
#' \eqn{\Lambda = (det(S) / trace((1/p)*S)^p})^{n/2},
#' where S is MLE of Sigma based on sample size n from \eqn{X ~ N_p(\mu,\Sigma)}.

#' Under null hypothesis, distribution of the test statistic \eqn{\Lambda} is
#' \eqn{\Lambda ~  prod_{j=2}^{p} (B_j)^{n/2}},
#' with \eqn{B_j ~ Beta((n-j)/2,(j-1)/p + (j-1)/2)}, where \eqn{j = (2,...,p)}.
#' Here we assume that \eqn{n > p}.

#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic \eqn{\Lambda}, say \eqn{W = -log(\Lambda)} is given by
#' \eqn{cf =  function(t) {cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j-1)/p + (j-1)/2)}},
#' where \eqn{j = (2,..., p)'}.
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
#' @note Ver.: 16-Sep-2018 21:11:16 (consistent with Matlab CharFunTool v1.3.0, 20-Jan-2018 12:43:15).
#'
#' @example R/Examples/example_LRT05_Sphericity.R
#'
#' @export
#'
LRT05_Sphericity <- function(W, n, p, options) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(n) || missing(p)) {
    stop("Enter input parameters n, p.")
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

  # CHARACTERISTIC FUNCTION CF
  coef <- options$coef
  if (length(coef) == 0) {
    coef <- -n / 2      # set this option for using with W = -log(LRT)
    # coef <- -1      # set this options for using normalized -log(LRT^(2/n))
  }

  ind <- Conj(2:p)
  alpha <- (n - ind) / 2
  beta <- (ind - 1) / p + (ind - 1) / 2
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
