#' @title
#' p-value of the log-transformed LRT statistic and/or its null
#' distribution CF/PDF/CDF
#'
#' @description
#' \code{LRT03_EqualityCovariances(W, n, p, q, options)} computes \eqn{p}-value of the log-transformed LRT statistic,
#' \eqn{W = -log(\Lambda)}, for testing the null hypothesis of equality
#' of covariance matrices (under normality assumptions) of \eqn{q} (\eqn{q > 1})
#' \eqn{p}-dimensional populations, and/or its null distribution CF/PDF/CDF.
#'
#' This is based on BALANCED samples of size n for each population!
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
#' In particular, let \eqn{X_k ~ N_p(mu_k,\Sigma_k}) are \eqn{p}-dimensional random
#' vectors, for \eqn{k = 1,...,q}. We want to test the hypothesis that the
#' covariance matrix \eqn{\Sigma} is common for all \eqn{X_k, k = 1,...,q}. Then, the
#' null hypothesis is given as \eqn{H0: \Sigma_1 = ... = \Sigma_q},
#' i.e. the covariance matrices are equal in all \eqn{q} populations. Here, the
#' LRT test statistic is given by
#' \eqn{\Lambda = ( q^{p*q} * prod(det(S_k)) / (det(S))^q )^{n/2}},
#' where \eqn{S_k} are MLEs of \eqn{\Sigma_k}, for \eqn{k = 1,...,q}, and \eqn{S = S_1 + ... + S_q},
#' based on \eqn{n} samples from each of the the \eqn{q} \eqn{p}-dimensional populations.
#'
#' Under null hypothesis, distribution of the test statistic \eqn{\Lambda} is
#' \eqn{\Lambda ~  prod_{k=1}^q prod_{j=1}^{p} (B_{jk})^{n/2}},
#' with \eqn{B_{jk} ~ Beta((n-j)/2,(j*(q-1)+2*k-1-q)/2)}, and we set \eqn{B_{11} = 1}
#' for \eqn{j=k=1}. Here we assume that  \eqn{n > p}.
#'
#' Hence, the exact characteristic function of the null distribution of
#' minus log-transformed LRT statistic \eqn{\Lambda}, say \eqn{W = -log(\Lambda)} is given by
#' \eqn{cf = function(t) {cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q))}},
#' where \eqn{k = (1*o,...,q*o)} with \eqn{p}-dimensional vector of ones \eqn{o = (1,...,1)}
#' and \eqn{j = (j_1,...,j_q)} with \eqn{j_k = 1:p}.
#'
#' @return
#' \eqn{p}-value of the log-transformed LRT statistic, \eqn{W = -log(\Lambda)}
#' and/or its null distribution CF/PDF/CDF.
#'
#' @references
#' [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate Statistical Analysis.
#' New York: Wiley, 3rd Ed., 2003. \cr
#' [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
#' near-exact distribution theory for the most common likelihood ratio
#' test statistics used in Multivariate Analysis. Test, 2011, 20.1:180-203. \cr
#' [3] WITKOVSKY, Viktor. Exact distribution of selected multivariate test
#' criteria by numerical inversion of their characteristic functions.
#' arXiv preprint arXiv:1801.02248, 2018.
#'
#' @family Likelihood Ratio Test
#'
#' @example R/Examples/example_LRT03_EqualityCovariances.R
#'
#' @export
#'
LRT03_EqualityCovariances <- function(W, n, p, q, options) {
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

  ##

  if (n <= p) {
    stop("Sample size n is too small.")
  }

  alpha <- rep(0, p * q)
  beta <- rep(0, p * q)
  ind <- 0
  for (k in 1:q) {
    for (j in 1:p) {
      ind <- ind + 1
      alpha[ind] <- (n - j) / 2
      beta[ind] <- (j * (q - 1) + 2 * k - 1 - q) / (2 * q)
    }
  }

  # For k=j=1 the coefficient beta=0, hence the term log(Beta(alpha,beta))=0.
  alpha <- alpha[2:length(alpha)]
  beta <- beta[2:length(beta)]

  # CHARACTERISTIC FUNCTION CF
  coef <- options$coef
  if (length(coef) == 0) {
    coef <- -n / 2    # set this option for using with W = -log(LRT)
    # coef <- -1    # set this options for using normalized -log(LRT^(2/n))
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
