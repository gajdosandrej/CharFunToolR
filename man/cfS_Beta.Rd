% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfS_Beta.R
\name{cfS_Beta}
\alias{cfS_Beta}
\title{Characteristic function of the zero-mean symmetric BETA distribution}
\usage{
cfS_Beta(t, theta = 1, coef, niid)
}
\arguments{
\item{t}{vector or array of real values, where the CF is evaluated.}

\item{theta}{the 'shape' parameter \code{theta > 0}. If empty, default value is \code{theta = 1}.}

\item{coef}{vector of coefficients of the linear combination of Beta distributed random variables.
If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.}

\item{niid}{scalar convolution coeficient.}
}
\value{
Characteristic function \eqn{cf(t)} of the Beta distribution.
}
\description{
\code{cfS_Beta(t, theta, coef, niid)} evaluates the characteristic function \eqn{cf(t)} of
the zero-mean symmetric BETA distribution defined on the interval \eqn{(-1,1)}.

\code{cfS_Beta} is an ALIAS of the more general function \code{cf_BetaSymmetric},
used to evaluate the characteristic function of a linear combination
of independent BETA distributed random variables.

The characteristic function of \eqn{X ~ BetaSymmetric(\theta)} is defined by
\deqn{cf(t) = cf_BetaSymmetric(t,\theta) = gamma(1/2+\theta) * (t/2)^(1/2-\theta) * besselj(\theta-1/2,t).}
}
\note{
Ver.: 16-Sep-2018 19:07:26 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
}
\examples{
## EXAMPLE1
# CF of the symmetric Beta distribution with theta = 3/2 on (-1,1)
theta <- 3 / 2
t <- seq(-50, 50, length.out = 501)
plotReIm(function(t)
        cfS_Beta(t, theta), t, title = "CF of the symmetric Beta distribution on (-1,1)")

## EXAMPLE2
# PDF/CDF of the the symmetric Beta distribution on (-1,1)
theta <- 3 / 2
cf <- function(t)
        cfS_Beta(t, theta)
x <- seq(-1, 1, length.out = 101)
xRange <- 2
options <- list()
options$dx <- 2 * pi / xRange
options$N <- 2 ^ 8
result <- cf2DistGP(cf = cf, x = x, options = options)
}
\references{
WITKOVSKY V. (2016). Numerical inversion of a characteristic
function: An alternative tool to form the probability distribution
of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Beta_distribution}.

Other Continuous Probability Distribution: 
\code{\link{cfS_Arcsine}()},
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
\code{\link{cfX_InverseGamma}()},
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

Other Symmetric Probability Distribution: 
\code{\link{cfS_Arcsine}()},
\code{\link{cfS_Gaussian}()},
\code{\link{cfS_Laplace}()},
\code{\link{cfS_Rectangular}()},
\code{\link{cfS_Student}()},
\code{\link{cfS_Trapezoidal}()},
\code{\link{cf_ArcsineSymmetric}()},
\code{\link{cf_BetaSymmetric}()},
\code{\link{cf_RectangularSymmetric}()},
\code{\link{cf_TSPSymmetric}()},
\code{\link{cf_TrapezoidalSymmetric}()}
}
\concept{Continuous Probability Distribution}
\concept{Symmetric Probability Distribution}
