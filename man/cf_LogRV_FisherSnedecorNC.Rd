% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_LogRV_FisherSnedecorNC.R
\name{cf_LogRV_FisherSnedecorNC}
\alias{cf_LogRV_FisherSnedecorNC}
\title{Characteristic function of a linear combinationof independent
LOG-TRANSFORMED non-central Fisher-Snedecor random variables}
\usage{
cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef, niid, tol)
}
\arguments{
\item{t}{vector or array of real values, where the CF is evaluated.}

\item{df1}{vector of the  degrees of freedom \code{df1 > 0}. If empty, default value is \code{df1 = 1}.}

\item{df2}{vector of the  degrees of freedom \code{df2 > 0}. If empty, default value is \code{df2 = 1}.}

\item{delta}{vector of non-centrality parameters.}

\item{coef}{vector of the coefficients of the linear combination of the Beta distributed random variables.
If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.}

\item{niid}{scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.}

\item{tol}{tolerance factor for selecting the Poisson weights, i.e. such that \eqn{PoissProb > tol}.
If empty, default value is \code{tol = 1e-12}.}
}
\value{
Characteristic function \eqn{cf(t)} of a linear combination
of independent LOG-TRANSFORMED non-central Fisher-Snedecor random variables.
}
\description{
\code{cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef, niid, tol)} evaluates characteristic function
of a linear combination (resp. convolution) of independent LOG-TRANSFORMED non-central Fisher-Snedecor
random variables, with distributions \eqn{F(df1_i,df2_i,\delta_i)}.

That is, \code{cf_LogRV_FisherSnedecorNC} evaluates the characteristic function
\eqn{cf(t)} of  \eqn{Y = coef_i*log(X_1) +...+ coef_N*log(X_N)}, where \eqn{X_i ~ F(df1_i,df2_i,\delta_i)}
are inedependent RVs, with \eqn{df1_i} and \eqn{df2_i} degrees of freedom, and the noncentrality parameters \eqn{delta_i >0},
for \eqn{i =  1,...,N}.

The characteristic function of \eqn{Y = log(X) with X ~ F(df1,df2,\delta)} is Poisson mixture
of the CFs of the shifted log-transformed central F RVs of the form
\deqn{cf(t) = cf_LogRV_FisherSnedecorNC(t,df1,df2,\delta) = exp(-\delta/2) sum_{j=1}^Inf (\delta/2)^j/j! *
exp(1i*t*(df1+2*j)/df1) * cf_LogRV_FisherSnedecor(t,df1+2*j,df2),}
where cf_LogRV_FisherSnedecor(t,df1,df2) denotes CF of log-transformed centrally distributed
F RVs with parameters df1 and df2. For more details on  the non-central
F distribution see \code{cf_FisherSnedecorNC}.
Alternatively,
\deqn{cf(t) = (df2/df1)^(1i*t) * gamma(df1/2 + 1i*t) / gamma(df1/2) * gamma(df2/2 - 1i*t) / gamma(df2/2) * 1F1(-1i*t;df1/2;-delta/2),}
where \eqn{1F1(a;b;z)} is the confluent hypergeometric function, also known as the Kummer function \eqn{M(a,b,z)}.
Hence,the characteristic function of \eqn{Y  = coef(1)*Y1 + ... + coef(N)*YN}
is  \eqn{cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t)}, where \eqn{cf_Yi(t)}
is evaluated with the parameters \eqn{df1_i}, \eqn{df2_i}, and \eqn{delta_i}.
}
\note{
Ver.: 20-Sep-2018 19:44:50 (consistent with Matlab CharFunTool v1.3.0, 10-Aug-2018 15:46:49).
}
\examples{
## EXAMPLE 1
# CF of the log-transformed non-central F RV with delta = 1 and coef = -1
df1 <- 3
df2 <- 5
delta <- 1
coef  <- -1
t <- seq(from = -10,
         to = 10,
         length.out = 201)
plotReIm(function(t)
        cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef),
        t,
        title = 'CF of minus log-transformed F RV')

## EXAMPLE 2
# CDF/PDF of the minus log-transformed non-central F RV with delta = 1
df1 <- 3
df2 <- 5
delta <- 1
coef <- -1
cf <-
        function(t)
                cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef)
options <- list()
options$N  <- 2 ^ 12
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3
# CDF/PDF of the linear combination of log-transformed non-central F RVs
df1 <- c(5, 4, 3)
df2 <- c(3, 4, 5)
delta <- c(0, 1, 2)
coef  <- -1 / 3
cf <-
        function(t)
                cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef)
options <- list()
options$N  <- 2 ^ 12
result <- cf2DistGP(cf = cf, options = options)
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Noncentral_F-distribution}.

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

Other Non-central Probability Distribution: 
\code{\link{cf_BetaNC}()},
\code{\link{cf_FisherSnedecorNC}()},
\code{\link{cf_LogRV_BetaNC}()},
\code{\link{cf_LogRV_ChiSquareNC}()}
}
\concept{Continuous Probability Distribution}
\concept{Non-central Probability Distribution}
