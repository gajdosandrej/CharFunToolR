% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfN_Binomial.R
\name{cfN_Binomial}
\alias{cfN_Binomial}
\title{Characteristic function of the Binomial distribution}
\usage{
cfN_Binomial(t, n = 10, p = 1/2, cfX)
}
\arguments{
\item{t}{vector or array of real values, where the CF is evaluated.}

\item{n}{number of trials.}

\item{p}{success probability, \eqn{0 \le} \code{p} \eqn{\le 1}, default value \code{p = 1/2}.}

\item{cfX}{function.}
}
\value{
Characteristic function \eqn{cf(t)} of the Binomial distribution.
}
\description{
\code{cfN_Binomial(t, n, p)} evaluates the characteristic function \eqn{cf(t)} of the
Binomial distribution with the parameters \code{n} (number of trials, \code{n} is a natural number)
and \code{p} (success probability, \code{p} in \eqn{[0,1]}), i.e.
\deqn{cfN_Binomial(t, n, p) = (1 - p + p*exp(1i*t))^n.}
For more details see [4].

\code{cfN_Binomial(t, n, p, cfX)} evaluates the compound characteristic function
\deqn{cf(t) = cfN_Binomial(-1i*log(cfX(t)), n, p),} where \code{cfX} is function
handle of the characteristic function \eqn{cfX(t)} of a continuous distribution
and/or random variable \eqn{X}.

Note that such CF is characteristic function of the compound distribution,
i.e. distribution of the random variable \eqn{Y = X_1 + ... + X_N}, where \eqn{X_i ~ F_X}
are i.i.d. random variables with common CF \eqn{cfX(t)}, and \eqn{N ~ F_N} is
independent RV with its CF given by \eqn{cfN(t)}.
}
\note{
Ver.: 16-Sep-2018 18:57:37 (consistent with Matlab CharFunTool v1.3.0, 15-Nov-2016 13:36:26).
}
\examples{
## EXAMPLE1
# CF of the Binomial distribution with n = 25, p = 0.3
n <- 25
p <- 0.3
t <- seq(-15, 15, length.out = 1001)
plotReIm(function(t)
        cfN_Binomial(t, n, p),
        t,
        title = "CF of the Binomial distribution with n = 25, p = 0.3")

## EXAMPLE2
# CF of the compound Binomial-Exponential distribution
n <- 25
p <- 0.3
lambda <- 10
cfX <- function(t)
        cfX_Exponential(t, lambda)
t <- seq(-10, 10, length.out = 501)
plotReIm(function(t)
        cfN_Binomial(t, n, p, cfX),
        t,
        title = "CF of the compound Binomial-Exponential distribution")

## EXAMPLE3
# PDF/CDF of the compound Binomial-Exponential distribution
n <- 25
p <- 0.3
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
x <- seq(0, 5, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
}
\references{
[1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
distribution based on numerical inversion of the compound empirical
characteristic function of frequency and severity. Preprint submitted
to Insurance: Mathematics and Economics.

[2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM
for computing distributions of collective risk models. Preprint submitted
to Journal of Statistical Software.

[3] WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
An alternative tool to form the probability distribution
of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.

[4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN 3-87773-025-6.
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Binomial_distribution}.

Other Discrete Probability Distribution: 
\code{\link{cfN_Delaporte}()},
\code{\link{cfN_GeneralizedPoisson}()},
\code{\link{cfN_Geometric}()},
\code{\link{cfN_Logarithmic}()},
\code{\link{cfN_NegativeBinomial}()},
\code{\link{cfN_Poisson}()},
\code{\link{cfN_PolyaEggenberger}()},
\code{\link{cfN_Quinkert}()},
\code{\link{cfN_Waring}()}
}
\concept{Discrete Probability Distribution}
