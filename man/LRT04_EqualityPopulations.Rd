% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRT04_EqualityPopulations.R
\name{LRT04_EqualityPopulations}
\alias{LRT04_EqualityPopulations}
\title{p-value of the log-transformed LRT statistic and/or its null
distribution CF/PDF/CDF}
\usage{
LRT04_EqualityPopulations(W, n, p, q, options)
}
\arguments{
\item{W}{observed value of the minus log-transformed LRT statistic
\eqn{W = -log(\Lambda)}. If empty, the  algorithm evaluates the
CF/PDF/CDF and the quantiles of the null distribution of \eqn{W}.}

\item{n}{sample size, \eqn{n > min(p+q-1)}.}

\item{p}{common dimension of the vectors \eqn{X_k, k = 1,...q}.}

\item{q}{number of normal populations, \eqn{q > 1}.}

\item{options}{option structure, for more details see \code{\link{cf2DistGP}}. Moreover, \cr
\code{x} set vector of values where PDF/CDF is evaluated, \cr
\code{prob} set vector of probabilities for the quantiles, \cr
\code{coef}  set arbitrary multiplicator of the argument \code{t}
of the characteristic function. If empty, default value is \eqn{-n/2}
(standard value for minus log-transform of LRT). Possible
alternative is e.g. \code{coef = -1}, leading to \eqn{W = -(2/n)*log(LRT)}.}
}
\value{
\eqn{p}-value of the log-transformed LRT statistic, \eqn{W = -log(\Lambda)}
and/or its null distribution CF/PDF/CDF.
}
\description{
\code{LRT04_EqualityPopulations(W, n, p, q, options)} computes \eqn{p}-value of the log-transformed LRT
statistic \eqn{W = -log(\Lambda)}, for testing the null hypothesis of equality
of \eqn{q} (\eqn{q > 1}) \eqn{p}-dimensional normal populations, and/or its null distribution CF/PDF/CDF.
This is based on BALANCED samples of size n for each population!

This is an EXPERIMENTAL version. Correctness should be checked again!
}
\details{
In particular, let \eqn{X_k ~ N_p(\mu_k,\Sigma_k)}, for \eqn{k = 1,...,q}. We want to
test the hypothesis that the \eqn{q} normal populations are equally
distributed. That is, we want to test that the
mean vectors \eqn{\mu_k} are equal for all \eqn{k = 1,...,q}, as well as the
covariance matrices \eqn{\Sigma_k} are equal for all \eqn{k = 1,...,q}. Then,
the null hypothesis is given as \eqn{H0: \mu_1 = ... = \mu_q} & \eqn{\Sigma_1 = ... = \Sigma_k}.
Here, the null hypothesis H0 and the LRT statistic can be decomposed:
\eqn{\Lambda = \Lambda_Means * \Lambda_Covariances}
where (first) \eqn{\Lambda_Covariances} represents the LRT for testing equality of
covariance matrices of given \eqn{q} normal populations, and (second)
\eqn{\Lambda_Means} represents (conditionally) the LRT for testing equality of
means of given \eqn{q} normal populations.

Under null hypothesis, distributions of \eqn{\Lambda_Covariances} and
\eqn{\Lambda_Means} are independent, and the distribution of the test statistic
\eqn{\Lambda} is \eqn{Lambda ~  \Lambda_Means * \Lambda_Covariances},
\eqn{~  (prod_{k=1}^q prod_{j=1}^{p} (B_{jk})^{n/2})* (prod_{j=1}^{p} (B_j)^{n*q/2})}
where the \eqn{B_{jk}} and \eqn{B_j} are mutually independent beta distributed
random variables. Here we assume that \eqn{n} is equal sample size for each
sample, \eqn{k = 1,...,q}, \eqn{n > p}.

Hence, the exact characteristic function of the null distribution of
minus log-transformed LRT statistic Lambda, say \eqn{W = -log(\Lambda)} is given by
\eqn{cf = function(t) {cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q))}}
. * cf_LogRV_Beta(-(n*q/2)*t, ((n-1)*q-i+1)/2, (q-1)/2),
where \eqn{i = (1, 2, ..., p)}, \eqn{k = (1*o,...,q*o)} with \eqn{p}-dimensional vector
of ones \eqn{o = (1,...,1)}  and \eqn{j = (j_1,...,j_q)} with \eqn{j_k = 1:p}.
}
\note{
Ver.: 16-Sep-2018 21:10:45 (consistent with Matlab CharFunTool v1.3.0, 20-Jan-2018 12:43:15).
}
\examples{
## EXAMPLE
# LRT for testing hypothesis on equality of populations
# Null distribution of the minus log-transformed LRT statistic
n <- 30                  # total sample size
p <- 8                   # dimension of X_k, k = 1,...,q where q = 5
q <- 5                   # number of populations
# W <- vector()          # observed value of W = -log(Lambda)
options <- list()
# options.coef = -1;
options$prob <- c(0.9, 0.95, 0.99)
output <- LRT04_EqualityPopulations(n = n, p = p, q = q, options = options)
str(output)
}
\references{
[1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate Statistical Analysis.
New York: Wiley, 3rd Ed., 2003.

[2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
near-exact distribution theory for the most common likelihood ratio
test statistics used in Multivariate Analysis. Test, 2011, 20.1:180-203.

[3] WITKOVSKY, Viktor. Exact distribution of selected multivariate test
criteria by numerical inversion of their characteristic functions. \cr
arXiv preprint arXiv:1801.02248, 2018.
}
\seealso{
Other Likelihood Ratio Test: 
\code{\link{LRT01_Independence}()},
\code{\link{LRT02_EqualityMeans}()},
\code{\link{LRT03_EqualityCovariances}()},
\code{\link{LRT05_Sphericity}()}
}
\concept{Likelihood Ratio Test}
