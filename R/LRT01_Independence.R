#' @title
#' p-value of the log-transformed LRT statistic and/or its null
#' distribution CF/PDF/CDF
#'
#' @description
#' \code{LRT01_Independence(W, n, p, options)} computes \eqn{p}-value of the log-transformed LRT statistic,
#' \eqn{W = -log(\Lambda)}, for testing the null hypothesis of independence (under
#' normality assumptions) of \eqn{m} groups of variables (\eqn{m > 1}), and/or its null
#' distribution CF/PDF/CDF.
#'
#' @param W observed value of the minus log-transformed LRT statistic,
#' \eqn{W = -log(\Lambda)}. If empty, the  algorithm evaluates the
#' CF/PDF/CDF and the quantiles of the null distribution of \eqn{W}.
#' @param n sample size (\eqn{n > q_k}).
#' @param p vector \eqn{p = (p_1,...,p_k)} of dimensions of \eqn{X_k, k = 1,...m}.
#' @param option option structure, for more details see \code{\link{cf2DistGP}}.
#' Moreover, \cr
#' \code{x} set vector of values where PDF/CDF is evaluated, \cr
#' \code{prob} set vector of probabilities for the quantiles, \cr
#' \code{coef} set arbitrary multiplicator of the argument \code{t}
#' of the characteristic function. If empty, default value is \eqn{-n/2}
#' (standard value for minus log-transform of LRT). Possible
#' alternative is e.g. \code{coef = -1}, leading to \eqn{W = -(2/n)*log(LRT)}.
#'
#' @details
#' In particular, let \eqn{X_k ~ N_{p_k}(\mu_k,\Sigma_k)}
#' are \eqn{p_k} dimensional random vectors, for \eqn{k = 1,...,m}. Let us denote
#' \eqn{X = (X_1,...,X_m)} and assume \eqn{X ~ N_p(mu,\Sigma)}. Then, the null hypothesis is
#' given as \eqn{H0: \Sigma = diag(\Sigma_1,...,\Sigma_m)},
#' i.e. the off-diagonal blocks of Sigma are blocks of zeros. Here, the LRT
#' test statistic is given by
#' \eqn{\Lambda = det(S) / prod(det(S_k))},
#' where \eqn{S} is MLE of Sigma, and \eqn{S_k} are MLEs of \eqn{\Sigma_k}, for \eqn{k = 1,...,m},
#' based on n samples from the compound vector \eqn{X = (X_1,...,X_m)}.
#'
#' Under null hypothesis, distribution of the test statistic \eqn{\Lambda} is
#' \eqn{\Lambda ~  prod_{k=1}^{m-1} prod_{j=1}^{p_k} (B_{jk})^{n/2}},
#' with \eqn{B_{jk} ~ Beta((n-q_k-j)/2,q_k/2)}, where \eqn{q_k = p_{k+1} + ... + p_m}.
#' Here we assume that \eqn{n > q_k}, for all \eqn{k = 1,...,m-1}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic \eqn{\Lambda}, say \eqn{W = -log(\Lambda)} is
#' given by \eqn{cf = function(t) {cf_LogRV_Beta(-(n/2)*t, (n-q-j)/2, q/2)}},
#' where \eqn{q = (q_1,...,q_m)} with \eqn{q_k = p_{k+1} + ... + p_m},
#' and \eqn{j = (j_1,...,j_m)} with \eqn{j_k = 1:p_k}.
#'
#' @return
#' \eqn{p}-value of the log-transformed LRT statistic, \eqn{W = -log(\Lambda)}
#' and/or its null distribution CF/PDF/CDF.
#'
#' @references
#' [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate Statistical Analysis.
#' New York: Wiley, 3rd Ed., 2003.\cr
#' [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
#' near-exact distribution theory for the most common likelihood ratio
#' test statistics used in Multivariate Analysis. Test, 2011, 20.1:180-203. \cr
#' [3] WITKOVSKY, Viktor. Exact distribution of selected multivariate test
#' criteria by numerical inversion of their characteristic functions. \cr
#' arXiv preprint arXiv:1801.02248, 2018.
#'
#' @family Likelihood Ratio Test
#'
#' @example R/Examples/example_LRT01_Independence.R
#'
#' @export
#'
LRT01_Independence <- function(W, n, p, options) {
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

  ##

  m <- length(p)
  q <- cumsum(p)
  N <- q[m - 1]

  if (n <= N) {
    stop("Sample size n is too small.")
  }

  alpha <- rep(0, N)
  beta <- rep(0, N)
  ind <- 0
  for (k in 1:(m - 1)) {
    for (j in 1:(p[k])) {
      ind <- ind + 1
      alpha[ind] <- (n - q[k] - j) / 2
      beta[ind] <- q[k] / 2
    }
  }

  # CHARACTERISTIC FUNCTION CF
  coef <- options$coef
  if (length(coef) == 0) {
    coef <- -n / 2  # set this option for using with W = -log(LRT)
    #coef <- -1   # set this options for using normalized -log(LRT^(2/n))
  }
  cf <- function(t) {
    cf_LogRV_Beta(t, alpha, beta, coef)
  }

  # Evaluate the p-value and PDF/CDF/QF of the log-transformed LRT statistic
  if (!missing(W) && length(W) > 0) {
    # P-VALUE
    options$isPlot <- FALSE
    result <- cf2DistGP(cf = cf, x = W, option = options)
    pval <- 1 - result$cdf
  } else {
    # DISTRIBUTION of Lambda PDF/CDF
    pval <- vector()
    result <-
      cf2DistGP(
        cf = cf,
        x = options$x,
        prob = options$prob,
        option = options
      )
  }

  # Save the parameters of the used beta distributions
  result$alpha <- alpha
  result$beta <- beta

  return(list(pval, result))
}

