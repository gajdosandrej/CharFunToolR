% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_TSPSymmetric.R
\name{cf_TSPSymmetric}
\alias{cf_TSPSymmetric}
\title{Characteristic function of a linear combination (resp. convolution) of
independent symmetric (location and scale shifted) TWO-SIDED-POWER (TSP)
random variables.}
\usage{
cf_TSPSymmetric(t, theta, mu, sigma, coef, niid)
}
\arguments{
\item{t}{vector or array of real values, where the CF is evaluated.}

\item{theta}{vector of the shape parameters \code{theta > 0}. If theta is scalar, it is assumed that all parameters theta are equal.
If empty, default value is \code{theta = 1}.}

\item{mu}{vector of location parameters, mu in Real. If empty, default value is \code{mu = 0}.}

\item{sigma}{vector of scale parameters, \code{sigma_i > 0}. If empty, default value is \code{sigma = 1}.}

\item{coef}{vector of the coefficients of the linear combination of the log-transformed random variables. If coef is scalar, it is
assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.}

\item{niid}{scalar convolution coeficient n, such that \eqn{Z = Y + ... + Y}
is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)} is independently
and identically distributed random variable. If empty, default value is \code{niid = 1}.}
}
\value{
Characteristic function of a linear combination (resp. convolution) of independent symmetric (location and scale shifted) TWO-SIDED-POWER (TSP) random variables.
}
\description{
That is, \code{cfS_TSPSymmetric(t, theta, mu, sigma, coef, niid)} evaluates the characteristic function \eqn{cf(t)} of
\eqn{Y = sum_{i=1}^N coef_i * (mu_i + sigma_i * X_i)}, where \eqn{X_i ~
TSP(theta_i)} are inedependent RVs, with symmetric TSP distributions
defined on the interval \eqn{(-1,1)} with zero mean and variance \eqn{Var(X_i) =
2*theta_i*gamma(theta_i)/gamma(3+theta_i)}, where \code{theta_i > 0} are shape
parameters for all zeqn{i = 1,...,N}.

 The characteristic function of the random variable\eqn{mu + sigma*X}, where
  X ~ TSP(theta) is given by
  \code{cf(t) = cfS_TSPSymmetric(t,theta,mu,sigma) = 1/2 * exp(1i*t*mu) * (hypergeom1F1(1,1+theta,1i*t*sigma) + hypergeom1F1(1,1+theta,-1i*t*sigma))}.

  Hence, the characteristic function of \eqn{Y  = coef_1*(mu_1+sigma_1*X_1) + coef_N*(mu_N+sigma_N*X_N)} is \eqn{cf_Y(t) = exp(1i*mu*t) *(cf_1(coef_1*sigma_1*t) * cf_N(coef_N*sigma_N*t))}, where \eqn{cf_i(t)} is
  the characteristic function of \eqn{X_i ~ TSP(theta_i)}.

SPECIAL CASES: \cr
1) \eqn{\theta = 1/2}; Arcsine distribution on \eqn{(-1,1)}:     \eqn{cf(t) = besselj(0,t)}, \cr
2) \eqn{\theta = 1};   Rectangular distribution on \eqn{(-1,1)}: \eqn{cf(t) = sin(t)/t}.
}
\note{
Ver.: 11-Aug-2021 16:18:40 (consistent with Matlab CharFunTool v1.5.1, 24-Jun-2017 18:25:56).
}
\examples{
## EXAMPLE 1
# CF of the symmetric TSP distribution with theta =3/2 on (-1,1)
theta <- 3/2
t <- seq(from = -50,
         to = 50,
         length.out =501)
plotReIm(function(t)
        cf_TSPSymmetric(t, theta),
        t,
        title = "CF of the symmetric TSP distribution  on (-1,1)")



##EXAMPLE2
# PDF/CDF of the symmetric TSP distribution on (-1,1)
thet <- 3/2
cf <- function(t)
        cf_TSPSymmetric(t, theta)
x <- seq(-1,1,length.out = 101)
xRange <- 2
options <- list()
options.N <-2^8
options.dt <- 2*pi/xRange
result <- cf2DistGP(cf, x, c(), options)

##EXAMPLE 3
# CF of the weighted linear combination of TSP RVs
theta <- c(1, 2, 3, 4, 5)/2
mu <- c(1, 2, 0, 0, 0)
sigma <- c(1, 2, 3, 4, 5)/5
coef <- 1/5
t <- seq(-50,50, length.out=  501)
plotReIm(function(t)
        cf_TSPSymmetric(t, theta, mu, sigma, coef),
        t,
        title = "CF of the weighted linear combination of TSP RVs")

## EXAMPLE 4
# CDF/PDF of the weighted linear combination of TSP RVs
thet <- c(1, 2, 3, 4, 5)/2
mu <- 0
sigma <- c(5, 4, 3, 2, 1)
coef <- 1/5
t <-seq(-50, 50, length.out = 501)
cf <- function(t)
        cf_TSPSymmetric(t, theta, mu, sigma, coef)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options.N = 2^12
result <- cf2DistGP(cf, c(), prob, options)
}
\references{
[1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
 function: An alternative tool to form the probability distribution of
 output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
 [2] VAN DORP, R.J., KOTZ, S. (2003). Generalizations of two-sided power
 distributions and their convolution. Communications in
 Statistics-Theory and Methods, 32(9), 1703-1723.
}
\seealso{
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
\code{\link{cf_TrapezoidalSymmetric}()},
\code{\link{cf_TriangularSymmetric}()},
\code{\link{cf_vonMises}()}

Other Symmetric Probability Distribution: 
\code{\link{cfS_Arcsine}()},
\code{\link{cfS_Beta}()},
\code{\link{cfS_Gaussian}()},
\code{\link{cfS_Laplace}()},
\code{\link{cfS_Rectangular}()},
\code{\link{cfS_Student}()},
\code{\link{cfS_Trapezoidal}()},
\code{\link{cf_ArcsineSymmetric}()},
\code{\link{cf_BetaSymmetric}()},
\code{\link{cf_RectangularSymmetric}()},
\code{\link{cf_TrapezoidalSymmetric}()}
}
\concept{Continuous Probability Distribution}
\concept{Symmetric Probability Distribution}
