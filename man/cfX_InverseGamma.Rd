% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfX_InverseGamma.R
\name{cfX_InverseGamma}
\alias{cfX_InverseGamma}
\title{Characteristic function of the INVERSE GAMMA distribution}
\usage{
cfX_InverseGamma(t, alpha, beta, coef, niid)
}
\arguments{
\item{t}{vector or array of real values, where the CF is evaluated.}

\item{alpha}{the shape parameter \code{alpha > 0}. If empty, default value is \code{alpha = 1}.}

\item{beta}{the rate (1/scale) parameter \code{beta > 0}. If empty, default value is \code{beta = 1}.}

\item{coef}{vector of coefficients of the linear combination of Inverse Gamma distributed random variables.
If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.}

\item{niid}{scalar convolution coeficient.}
}
\value{
Characteristic function \eqn{cf(t)} of the INVERSE GAMMA distribution.
}
\description{
\code{cfX_InverseGamma(t, alpha, beta, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
of the INVERSE GAMMA distribution with the shape parameter \code{alpha > 0} and the rate parameter \code{beta > 0}.

\code{cfX_InverseGamma} is an ALIAS NAME of the more general function \code{cf_InverseGamma},
used to evaluate the characteristic function of a linear combination
of independent INVERSE GAMMA distributed random variables.

The characteristic function of the GAMMA distribution is defined
by \deqn{cf(t) =  2 / gamma(\alpha) * (-1i*\beta*t).^(\alpha/2) * besselk(\alpha,sqrt(-4i*\beta*t)).}
}
\note{
Ver.: 16-Sep-2018 19:25:22 (consistent with Matlab CharFunTool v1.3.0, 15-Nov-2016 13:36:26).
}
\examples{
## EXAMPLE 1
# CF of the InverseGamma distribution with alpha = 2, beta = 2
alpha <- 2
beta <- 2
t <- seq(-20, 20, length.out = 501)
plotReIm(function(t)
        cfX_InverseGamma(t, alpha, beta), t,
        title = "CF of the InverseGamma distribution with alpha = 2, beta = 2")

## EXAMPLE 2
# PDF/CDF of the InverseGamma distribution with alpha = 2, beta = 2
alpha <- 2
beta <- 2
cf <- function(t)
        cfX_InverseGamma(t, alpha, beta)
x <- seq(0, 15, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
options$SixSigmaRule <- 10
options$N <- 2 ^ 14
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE 3
# PDF/CDF of the compound Binomial-InverseGamma distribution
p <- 0.3
n <- 25
alpha <- 2
beta <- 2
cfX <- function(t)
        cfX_InverseGamma(t, alpha, beta)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
x <- seq(0, 70, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
}
\references{
WITKOVSKY, V.: Computing the distribution of a linear combination
of inverted gamma variables, Kybernetika 37 (2001), 79-90.
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}.

Other Continuous Probability Distribution: 
\code{\link{cfS_Arcsine}()},
\code{\link{cfS_Beta}()},
\code{\link{cfS_Gaussian}()},
\code{\link{cfS_Laplace}()},
\code{\link{cfS_Rectangular}()},
\code{\link{cfS_Student}()},
\code{\link{cfS_TSP}()},
\code{\link{cfS_Trapezoidal}()},
\code{\link{cfS_Triangular}()},
\code{\link{cfS_Wigner}()},
\code{\link{cfX_ChiSquare}()},
\code{\link{cfX_Exponential}()},
\code{\link{cfX_FisherSnedecor}()},
\code{\link{cfX_Gamma}()},
\code{\link{cfX_LogNormal}()},
\code{\link{cf_ArcsineSymmetric}()},
\code{\link{cf_BetaNC}()},
\code{\link{cf_BetaSymmetric}()},
\code{\link{cf_Beta}()},
\code{\link{cf_ChiSquare}()},
\code{\link{cf_Exponential}()},
\code{\link{cf_FisherSnedecorNC}()},
\code{\link{cf_FisherSnedecor}()},
\code{\link{cf_Gamma}()},
\code{\link{cf_InverseGamma}()},
\code{\link{cf_Laplace}()},
\code{\link{cf_LogRV_BetaNC}()},
\code{\link{cf_LogRV_Beta}()},
\code{\link{cf_LogRV_ChiSquareNC}()},
\code{\link{cf_LogRV_ChiSquare}()},
\code{\link{cf_LogRV_FisherSnedecorNC}()},
\code{\link{cf_LogRV_FisherSnedecor}()},
\code{\link{cf_LogRV_MeansRatioW}()},
\code{\link{cf_LogRV_MeansRatio}()},
\code{\link{cf_LogRV_WilksLambdaNC}()},
\code{\link{cf_LogRV_WilksLambda}()},
\code{\link{cf_Normal}()},
\code{\link{cf_RectangularSymmetric}()},
\code{\link{cf_Student}()},
\code{\link{cf_TSPSymmetric}()},
\code{\link{cf_TrapezoidalSymmetric}()},
\code{\link{cf_TriangularSymmetric}()},
\code{\link{cf_vonMises}()}
}
\concept{Continuous Probability Distribution}
