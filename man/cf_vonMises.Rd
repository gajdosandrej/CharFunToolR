% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_vonMises.R
\name{cf_vonMises}
\alias{cf_vonMises}
\title{Characteristic function of a linear combination
of independent VON MISES random variables}
\usage{
cf_vonMises(t, mu = 0, kappa = 1, coef, niid)
}
\arguments{
\item{t}{numerical values (number, vector...).}

\item{mu}{in \eqn{(-\pi, \pi)}.}

\item{kappa}{\eqn{> 0}.}

\item{coef}{vector of the coefficients of the linear combination
of the IGamma random variables. If coef is scalar, it is assumed
that all coefficients are equal. If empty, default value is \code{coef = 1}.}

\item{niid}{scalar convolution coeficient \code{niid}, such that \eqn{Z = Y +...+ Y} is
sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef_i * X_i}
is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.}
}
\value{
Characteristic function \eqn{cf(t)} of a linear combination
of independent VON MISES random variables.
}
\description{
\code{cf_vonMises(t, mu, kappa, coef, niid)} evaluates the characteristic function \eqn{cf(t)}
of \eqn{Y = sum_{i=1}^N coef_i * X_i} where \eqn{X_i ~ vonMises(\mu_i,\kappa_i)}
inedependent RVs, with the locarion parameters \eqn{\mu_i} in Real and the rate
parameters \eqn{\kappa_i > 0}, for \eqn{i = 1,...,N}.

The characteristic function of the \eqn{vonMises(\mu,\kappa)} distribution is
\deqn{cf(t) = cf_vonMises(t,\mu,\kappa) = besseli(t,\kappa)/besseli(0,\kappa) * exp(1i*t*\mu).}
}
\details{
The VON MISES distribution is circular distribution on the interval
of length \eqn{2*\pi}, here we consider \eqn{(-\pi,\pi)}, equivalent of the normal
distribution with the real parameter \eqn{\mu} and rate parameter \eqn{\kappa > 0}
(\eqn{\mu} and \eqn{1/\kappa} are analogous to \eqn{\mu} and \eqn{\sigma^2}, the mean and variance
in the normal distribution), on a whole circle, i.e. the interval of angles \eqn{(-\pi,\pi)}.
}
\note{
Ver.: 16-Sep-2018 18:41:05 (consistent with Matlab CharFunTool v1.3.0, 24-Jun-2017 18:25:56).
}
\examples{
## EXAMPLE 1:
# CF of the weighted linear combinantion of the von Mises RVs
mu <- c(0, 0,-pi / 2, pi / 2, 0)
kappa <- c(1, 2, 3, 4, 5)
coef  <- c(1, 2, 3, 4, 5) / 15
t <- seq(from = -20,
         to = 20,
         length.out = 201)
plotReIm(function(t)
        cf_vonMises(t,
                    mu,
                    kappa,
                    coef),
        t,
        title = 'CF of the weighted linear combinantion of the von Mises RVs')

## EXAMPLE2
# CDR/PDF of the weighted linear combinantion of the von Mises RVs
mu <- c(0, 0,-pi / 2, pi / 2, 0)
kappa <- c(1, 2, 3, 4, 5)
coef  <- c(1, 2, 3, 4, 5) / 15
t <- seq(from = -20,
         to = 20,
         length.out = 201)
cf <- function(t) {
        cf_vonMises(t, mu, kappa, coef)
}
x <- seq(from = -pi,
         to = pi,
         length.out = 201)
prob  <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- -pi
options$xMax <- pi
result <- cf2DistGP(cf, x, prob, options)
result
angle  <- result$x
radius <- result$pdf
plotPolar(angle, radius)
# figure; polarplot(angle,radius); ???
# x = gca; ax.ThetaAxisUnits = 'radians'; ???

## EXAMPLE3
# CF of the mixture of the von Mises distribution on (-pi,pi)
mu1 <- 0
kappa1 <- 5
mu2 <- 1
kappa2 <- 15
cf <-
        function(t) {
                0.25 * cf_vonMises(t, mu1, kappa1) + 0.75 * cf_vonMises(t, mu2, kappa2)
        }
options <- list()
options$xMin <- -pi
options$xMax <- pi
result <- cf2DistGP(cf = cf, options = options)
angle  <- result$x
radius <- result$pdf
plotPolar(angle, radius)
#   figure; polarplot(angle,radius);
#   ax = gca; ax.ThetaAxisUnits = 'radians';

## EXAMPLE4
# PDF/CDF of the mixture of the von Mises distribution on (0,2*pi)
mu1 <- 0
kappa1 <- 5
mu2 <- 1
kappa2 <- 15
mu3 <- pi
kappa3 <- 10
cf  <-
        function(t) {
                0.25 * cf_vonMises(t, mu1, kappa1) + 0.25 * cf_vonMises(t, mu2, kappa2) + 0.5 *
                        cf_vonMises(t, mu3, kappa3)
        }
options <- list()
options$isCircular   <- TRUE
options$correctedCDF <- TRUE
options$xMin <- 0
options$xMax <- 2 * pi
x <- seq(from = 0,
         to = 2 * pi,
         length.out = 100)
result <- cf2DistGP(cf = cf, x = x, options = options)
angle  <- result$x
radius <- result$pdf
plotPolar(angle, radius)
#figure; polarplot(angle,radius);
#ax = gca; ax.ThetaAxisUnits = 'radians';
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Von_Mises_distribution}.

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
\code{\link{cf_TriangularSymmetric}()}
}
\concept{Circular Probability Distribution}
\concept{Continuous Probability Distribution}
