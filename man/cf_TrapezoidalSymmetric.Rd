% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_TrapezoidalSymmetric.R
\name{cf_TrapezoidalSymmetric}
\alias{cf_TrapezoidalSymmetric}
\title{Characteristic function of a linear combination
of independent zero-mean symmetric TRAPEZOIDAL random variables}
\usage{
cf_TrapezoidalSymmetric(t, lambda, coef, niid)
}
\arguments{
\item{t}{vector or array of real values, where the CF is evaluated.}

\item{lambda}{parameter of the offset, \eqn{0 \le} \code{lambda} \eqn{\le 1}.
If empty, default value is \code{lambda = 0}}

\item{coef}{vector of the coefficients of the linear combination
of the zero-mean symmetric TRAPEZOIDAL random variables. If coef is scalar,
it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.}

\item{niid}{scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.}
}
\value{
Characteristic function \eqn{cf(t)} of a linear combination
of independent zero-mean symmetric TRAPEZOIDAL random variables.
}
\description{
\code{cf_TrapezoidalSymmetric(t, lambda, coef, niid)} evaluates
the characteristic function of a linear combination (resp. convolution)
of independent zero-mean symmetric TRAPEZOIDAL random variables defined on the interval \eqn{(-1,1)}.

That is, \code{cf_TrapezoidalSymmetric} evaluates the characteristic function
\eqn{cf(t)} of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ TrapezoidalSymmetric(\lambda_i)}
are independent RVs defined on \eqn{(-1,1)}, for all \eqn{i = 1,...,N}.

The characteristic function of \eqn{X ~ TrapezoidalSymmetric(\lambda)}, with
\eqn{E(X) = 0} and \eqn{Var(X) = (1+\lambda^2)/6}, is
\deqn{cf(t) = cf_TrapezoidalSymmetric(t) = cf_RectangularSymmetric(w*t))*cf_RectangularSymmetric((1-w)*t) = (sin(w*t)/(w*t))*(sin((1-w)*t)/((1-w)*t)),}
where \eqn{w  = (1+\lambda)/2}.
}
\note{
Ver.: 16-Sep-2018 18:39:29 (consistent with Matlab CharFunTool v1.3.0, 02-Jun-2017 12:08:24).
}
\examples{
## EXAMPLE 1
#CF of the symmetric Trapezoidal distribution with lambda = 1/2
lambda <- 1 / 2
t <- seq(from = -50,
         to = 50,
         length.out = 201)
plotReIm(function(t)
        cf_TrapezoidalSymmetric(t, lambda), t,
        title = "CF of the symmetric Trapezoidal distribution on (-1,1)")

## EXAMPLE 2
# CF of a linear combination of independent Trapezoidal RVs
t <- seq(from = -20,
         to = 20,
         length.out = 201)
lambda <- c(3, 3, 4, 4, 5) / 7
coef <- c(1, 2, 3, 4, 5) / 15
plotReIm(function(t)
        cf_TrapezoidalSymmetric(t, lambda, coef), t,
        title = "CF of a linear combination of independent Trapezoidal RVs")

## EXAMPLE 3
# PDF/CDF of a weighted linear combination of independent Trapezoidal RVs
lambda <- c(3, 3, 4, 4, 5) / 7
coef <- c(1, 2, 3, 4, 5) / 15
cf <- function(t)
        cf_TrapezoidalSymmetric(t, lambda, coef)
x <- seq(from = -1,
         to = 1,
         length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2 ^ 12
options$xMin <- -1
options$xMax <- 1
result <- cf2DistGP(cf, x, prob, options)
}
\references{
WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
An alternative tool to form the probability distribution
of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Trapezoidal_distribution}.

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
\code{\link{cf_TSPSymmetric}()}
}
\concept{Continuous Probability Distribution}
\concept{Symmetric Probability Distribution}
