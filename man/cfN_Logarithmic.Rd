% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfN_Logarithmic.R
\name{cfN_Logarithmic}
\alias{cfN_Logarithmic}
\title{Characteristic function of the Logarithmic distribution}
\usage{
cfN_Logarithmic(t, p = 1, cfX)
}
\arguments{
\item{t}{numerical values (number, vector...).}

\item{p}{success probability, \eqn{0 \le} \code{p} \eqn{\le 1}, default value \code{p = 1}.}

\item{cfX}{function.}
}
\value{
Characteristic function \eqn{cf(t)} of the Logarithmic distribution.
}
\description{
\code{cfN_Logarithmic(t, p, cfX)} evaluates the characteristic function \eqn{cf(t)} of the
Logarithmic distribution defined on non-negative integers \eqn{n =0, 1, \ldots},
with the parameters \code{p} (success probability \code{p} in \eqn{[0,1]}), i.e.
\deqn{cfN_Logarithmic(t, p) = log(1 - p * exp(1i*t)) / log(1 - p).}
For more details see [4].

\code{cfN_Logarithmic(t,p,cfX)} evaluates the compound characteristic function
\deqn{cf(t) = cfN_Logarithmic(-1i*log(cfX(t)), p),} where \code{cfX} is function
handle of the characteristic function \eqn{cfX(t)} of a continuous distribution
and/or random variable \eqn{X}.

Note that such CF is characteristic function of the compound distribution,
i.e. distribution of the random variable \eqn{Y = X_1 + ... + X_N}, where \eqn{X_i ~ F_X}
are i.i.d. random variables with common CF \eqn{cfX(t)}, and \eqn{N ~ F_N} is
independent RV with its CF given by \eqn{cfN(t)}.
}
\note{
Ver.: 16-Sep-2018 19:01:18 (consistent with Matlab CharFunTool v1.3.0, 15-Nov-2016 13:36:26).
}
\examples{
## EXAMPLE1
# CF of the Logarithmic distribution with the parameter p = 0.5
p <- 0.5
t <- seq(-10, 10, length.out = 501)
plotReIm(function(t)
        cfN_Logarithmic(t, p),
        t,
        title = "CF of the Logarithmic distribution with the parameter p = 0.5")

## EXAMPLE2
# CF of the compound Logarithmic-Exponential distribution
p <- 0.5
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
t <- seq(-10, 10, length.out = 501)
plotReIm(function(t)
        cfN_Logarithmic(t, p, cfX),
        t,
        title = "CF of the compound Logarithmic-Exponential distribution")

## EXAMPLE3
# PDF/CDF of the compound Logarithmic-Exponential distribution
p <- 0.5
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
cf <- function(t)
        cfN_Logarithmic(t, p, cfX)
x <- seq(0, 3, length.out = 101)
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
\url{https://en.wikipedia.org/wiki/Logarithmic_distribution}.

Other Discrete Probability Distribution: 
\code{\link{cfN_Binomial}()},
\code{\link{cfN_Delaporte}()},
\code{\link{cfN_GeneralizedPoisson}()},
\code{\link{cfN_Geometric}()},
\code{\link{cfN_NegativeBinomial}()},
\code{\link{cfN_Poisson}()},
\code{\link{cfN_PolyaEggenberger}()},
\code{\link{cfN_Quinkert}()},
\code{\link{cfN_Waring}()}
}
\concept{Discrete Probability Distribution}
