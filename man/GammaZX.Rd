% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GammaZX.R
\name{GammaZX}
\alias{GammaZX}
\title{Gamma function valid in the entire complex plane}
\usage{
GammaZX(z, funmode)
}
\arguments{
\item{z}{complex argument.}

\item{funmode}{function mode, \code{funmode = 0} for \eqn{ln[gamma(z)]} and \code{funmode = 1} for \eqn{gamma(z)}.}
}
\value{
Function returns \eqn{ln[gamma(z)]} or \eqn{gamma(z)}.
}
\description{
\code{GammaZX(z, funmode)} evaluates Gamma function valid in the entire complex plane, the argument \code{z}
may be complex and of any size.
}
\note{
Ver.: 01-Oct-2018 13:00:11 (consistent with Matlab CharFunTool v1.3.0, 20-Jul-2018 16:13:31).
}
\seealso{
Other Utility Function: 
\code{\link{ChebCoefficients}()},
\code{\link{ChebPoints}()},
\code{\link{ChebPolyValues}()},
\code{\link{ChebPoly}()},
\code{\link{ChebValues}()},
\code{\link{GammaLog}()},
\code{\link{GammaMultiLog}()},
\code{\link{GammaMulti}()},
\code{\link{Hypergeom1F1MatApprox}()},
\code{\link{Hypergeom1F1Mat}()},
\code{\link{Hypergeom2F1Mat}()},
\code{\link{Hypergeom2F1}()},
\code{\link{HypergeompFqMat}()},
\code{\link{InterpChebValues}()},
\code{\link{hypergeom1F1}()},
\code{\link{interpBarycentric}()}
}
\concept{Utility Function}
