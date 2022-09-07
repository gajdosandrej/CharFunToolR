#' @title
#' Evaluating CDF/PDF/QF from CF of a continous distribution F by using the Gil-Pelaez inversion formulae.
#'
#' @description
#' \code{cf2DistGP(cf, x, prob, options)} calcuates the CDF/PDF/QF from the Characteristic Function CF
#' by using the Gil-Pelaez inversion formulae:
#' \deqn{cdf(x) = 1/2 + (1/\pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt,}
#' \deqn{pdf(x) = (1/\pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.}
#'
#' The FOURIER INTEGRALs are calculated by using the simple TRAPEZOIDAL QUADRATURE method, see below.
#'
#' @family CF Inversion Algorithm
#'
#' @importFrom stats runif
#' @importFrom graphics plot grid
#'
#' @seealso For more details see:
#' \url{https://arxiv.org/pdf/1701.08299.pdf}.
#'
#' @param cf      function handle for the characteristic function CF.
#' @param x       vector of values where the CDF/PDF is computed.
#' @param prob    vector of values from \eqn{[0,1]} for which the quantile function is evaluated.
#' @param options  structure (list) with the following default parameters:
#' \itemize{
#'     \item \code{options$isCompound = FALSE} treat the compound distributions, of the RV \eqn{Y = X_1 + ... + X_N},
#'     where \eqn{N} is discrete RV and \eqn{X\ge0} are iid RVs from nonnegative continuous distribution,
#'     \item \code{options$isCircular = FALSE} treat the circular distributions on \eqn{(-\pi, \pi)},
#'     \item \code{options$isInterp = FALSE} create and use the interpolant functions for PDF/CDF/RND,
#'     \item \code{options$N = 2^10} N points used by FFT,
#'     \item \code{options$xMin = -Inf} set the lower limit of \code{x},
#'     \item \code{options$xMax = Inf} set the upper limit of \code{x},
#'     \item \code{options$xMean = NULL} set the MEAN value of \code{x},
#'     \item \code{options$xStd = NULL} set the STD value of \code{x},
#'     \item \code{options$dt = NULL} set grid step \eqn{dt = 2*\pi/xRange},
#'     \item \code{options$T = NULL} set upper limit of \eqn{(0,T)}, \eqn{T = N*dt},
#'     \item \code{options$SixSigmaRule = 6} set the rule for computing domain,
#'     \item \code{options$tolDiff = 1e-4} tol for numerical differentiation,
#'     \item \code{options$isPlot = TRUE} plot the graphs of PDF/CDF,
#'
#'     \item options$DIST                   list with information for future evaluations,
#'                                         \code{options$DIST} is created automatically after first call:
#'     \itemize{
#'         \item \code{options$DIST$xMin} the lower limit of \code{x},
#'         \item \code{options$DIST$xMax} the upper limit of \code{x},
#'         \item \code{options$DIST$xMean} the MEAN value of \code{x},
#'         \item \code{options$DIST$cft} CF evaluated at \eqn{t_j} : \eqn{cf(t_j)}.
#'         }
#'     }
#'
#' @return
#' \item{result}{structure (list) with CDF/PDF/QF amd further details.}
#'
#' @details
#' If \code{options.DIST} is provided, then \code{cf2DistGP} evaluates CDF/PDF based on
#' this information only (it is useful, e.g., for subsequent evaluation of
#' the quantiles). \code{options.DIST} is created automatically after first call.
#' This is supposed to be much faster, bacause there is no need for further
#' evaluations of the characteristic function. In fact, in such case the
#' function handle of the CF is not required, i.e. in such case do not specify \code{cf}.
#'
#' The required integrals are evaluated approximately by using the simple
#' trapezoidal rule on the interval \eqn{(0,T)}, where \eqn{T = N * dt} is a sufficienly
#' large integration upper limit in the frequency domain.
#'
#' If the optimum values of \eqn{N} and \eqn{T} are unknown, we suggest, as a simple
#' rule of thumb, to start with the application of the six-sigma-rule for
#' determining the value of \eqn{dt = (2*\pi)/(xMax-xMin)}, where \eqn{xMax = xMean +}
#' \eqn{6*xStd}, and \eqn{xMin = xMean - 6*xStd}, see \code{[1]}.
#'
#' Please note that THIS (TRAPEZOIDAL) METHOD IS AN APPROXIMATE METHOD:
#' Frequently, with relatively low numerical precision of the results of the
#' calculated PDF/CDF/QF, but frequently more efficient and more precise
#' than comparable Monte Carlo methods.
#'
#' However, the numerical error (truncation error and/or the integration
#' error) could be and should be properly controled!
#'
#' CONTROLING THE PRECISION:
#' Simple criterion for controling numerical precision is as follows: Set \eqn{N}
#' and \eqn{T = N*dt} such that the value of the integrand function
#' \eqn{Imag(e^(-1i*t*x) * cf(t)/t)} is sufficiently small for all \eqn{t > T}, i.e.
#' \eqn{PrecisionCrit = abs(cf(t)/t) <= tol},
#' for pre-selected small tolerance value, say \eqn{tol = 10^-8}. If this
#' criterion is not valid, the numerical precission of the result is
#' violated, and the method should be improved (e.g. by selecting larger \eqn{N}
#' or considering other more sofisticated algorithm - not considered here).
#'
#' @references
#' [1] WITKOVSKY, V.: On the exact computation of the density and
#' of the quantiles of linear combinations of t and F random variables.
#' Journal of Statistical Planning and Inference, 2001, 94, 1-13.
#'
#' [2] WITKOVSKY, V.: Matlab algorithm TDIST: The distribution
#' of a linear combination of Student's t random variables.
#' In COMPSTAT 2004 Symposium (2004), J. Antoch, Ed., Physica-Verlag/Springer 2004,
#' Heidelberg, Germany, pp. 1995-2002.
#'
#' [3] WITKOVSKY, V., WIMMER, G., DUBY, T. Logarithmic Lambert W x F
#' random variables for the family of chi-squared distributions
#' and their applications. Statistics & Probability Letters, 2015, 96, 223-231.
#'
#' [4] WITKOVSKY, V.: Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 2016, 5(3), 32-44.
#'
#' [5] WITKOVSKY, V., WIMMER, G., DUBY, T. Computing the aggregate loss distribution
#' based on numerical inversion of the compound empirical characteristic function
#' of frequency and severity. ArXiv preprint, 2017, arXiv:1701.08299.
#'
#' @note Ver.: 16-Sep-2018 17:59:07 (consistent with Matlab CharFunTool v1.3.0, 22-Sep-2017 11:11:11).
#'
#' @example R/Examples/example_cf2DistGP.R
#'
#' @export
#'
cf2DistGP <- function(cf, x, prob, options) {
    ## CHECK THE INPUT PARAMETERS
    start_time <- Sys.time()
    if (missing(x)) {
      x <- vector()
    }

    if (missing(prob)) {
      prob <- vector()
    }

    if (missing(options)) {
      options <- list()
    }

    if (is.null(options$isCompound)) {
      options$isCompound <- FALSE
    }

    if (is.null(options$isCircular)) {
      options$isCircular <- FALSE
    }

    if (is.null(options$N)) {
      if (options$isCompound) {
        options$N <- 2 ^ 14
      }  else {
        options$N <- 2 ^ 10
      }
    }

    if (is.null(options$xMin)) {
      if (options$isCompound) {
        options$xMin = 0
      }  else {
        options$xMin = -Inf
      }
    }

    if (is.null(options$xMax)) {
      options$xMax <- Inf
    }

    if (is.null(options$xMean)) {
      options$xMean <- vector()
    }

    if (is.null(options$xStd)) {
      options$xStd <- vector()
    }

    if (is.null(options$dt)) {
      options$dt <- vector()
    }

    if (is.null(options$T)) {
      options$T <- vector()
    }

    if (is.null(options$SixSigmaRule)) {
      if (options$isCompound) {
        options$SixSigmaRule <- 10
      }  else {
        options$SixSigmaRule <- 6
      }
    }

    if (is.null(options$tolDiff)) {
      options$tolDiff = 1e-04
    }

    if (is.null(options$crit)) {
      options$crit = 1e-12
    }

    if (is.null(options$isPlot)) {
      options$isPlot <- TRUE
    }

    if(is.null(options$DIST)) {
      options$DIST <- list()
    }

    # Other options parameters
    if (is.null(options$qf0)) {
      options$qf0 <- Re((cf(1e-4) - cf(-1e-4)) / (2e-4 * 1i))
    }

    if (is.null(options$maxiter)) {
      options$maxiter <- 1000
    }

    if (is.null(options$xN)) {
      options$xN <- 101
    }

    if (is.null(options$chebyPts)) {
      options$chebyPts <- 2 ^ 9
    }

    if (is.null(options$correctedCDF)) {
      if (options$isCircular) {
        options$correctedCDF <- TRUE
      } else {
        options$correctedCDF <- FALSE
      }
    }

    if (is.null(options$isInterp)) {
      options$isInterp <- FALSE
    }

    ## GET/SET the DEFAULT parameters and the OPTIONS
    cfOld <- vector()

    if (length(options$DIST) > 0) {
      xMean <- options$DIST$xMean
      cft <- options$DIST$cft
      xMin <- options$DIST$xMin
      xMax <- options$DIST$xMax
      cfLimit <- options$DIST$cfLimit
      xRange <- xMax - xMin
      dt <- 2 * pi / xRange
      N <- length(cft)
      t <- (1:N) * dt
      xStd <- numeric()
    } else {
      N <- options$N
      dt <- options$dt
      T <- options$T
      xMin <- options$xMin
      xMax <- options$xMax
      xMean <- options$xMean
      xStd <- options$xStd
      SixSigmaRule <- options$SixSigmaRule
      tolDiff <- options$tolDiff

      # Special treatment for compound distributions. If the real value of CF at infinity (large value)
      # is positive, i.e. cfLimit = real(cf(Inf)) > 0. In this case the
      # compound CDF has jump at 0 of size equal to this value, i.e. cdf(0) =
      # cfLimit, and pdf(0) = Inf. In order to simplify the calculations, here we
      # calculate PDF and CDF of a distribution given by transformed CF, i.e.
      # cf_new(t) = (cf(t)-cfLimit) / (1-cfLimit); which is converging to 0 at Inf,
      # i.e. cf_new(Inf) = 0. Using the transformed CF requires subsequent
      # recalculation of the computed CDF and PDF, in order to get the originaly
      # required values: Set pdf_original(0) =  Inf & pdf_original(x) = pdf_new(x) * (1-cfLimit),
      # for x > 0. Set cdf_original(x) =  cfLimit + cdf_new(x) * (1-cfLimit).

      cfLimit <- Re(cf(1e300))
      cfOld <- cf
      cf2 <- cf
      if (options$isCompound) {
        if (cfLimit > 1e-13) {
          # cf <- function(t) {
          #   (cf(t) - cfLimit) / (1 - cfLimit)
          cf2 <- function(t) {
            (cf(t) - cfLimit) / (1 - cfLimit)
          }
        }
        options$isNonnegative <- TRUE
      }

      cft <- cf2(tolDiff * (1:4))
      cftRe <- Re(cft)
      cftIm <- Im(cft)

      if (length(xMean) == 0) {
        if (options$isCircular) {
          xMean <- Arg(cf2(1))
        } else {
          xMean <-
            (8 * cftIm[1] / 5 - 2 * cftIm[2] / 5 + 8 * cftIm[3] / 105 - 2 * cftIm[4] /
               280) / tolDiff
        }
      }
      if (length(xStd) == 0) {
        if (options$isCircular) {
          # see https://en.wikipedia.org/wiki/Directional_statistics
          xStd <- sqrt(-2 * log(abs(cf2(1))))
        } else {
          xM2   <-
            (205 / 72 - 16 * cftRe[1] / 5 + 2 * cftRe[2] / 5 - 16 * cftRe[3] / 315 + 2 *
               cftRe[4] / 560) / (tolDiff ^ 2)
          xStd  <- sqrt(xM2 - xMean ^ 2)
        }
      }
      if (is.finite(xMin) && is.finite(xMax)) {
        xRange <- xMax - xMin
      } else if (length(dt) > 0) {
        xRange <- 2 * pi / dt
      } else if (length(T) > 0) {
        xRange <- 2 * pi / (T / N)
      } else {
        if (options$isCircular) {
          xMin <- -pi
          xMax <- pi
        } else {
          if (is.finite(xMin)) {
            xMax <- xMean + SixSigmaRule * xStd
          } else if (is.finite(xMax)) {
            xMin <- xMean - SixSigmaRule * xStd
          } else {
            xMin <- xMean - SixSigmaRule * xStd
            xMax <- xMean + SixSigmaRule * xStd
          }
        }
        xRange <- xMax - xMin
      }

      dt      <- 2 * pi / xRange
      t       <- (1:N) * dt
      cft     <- cf2(t)
      cft[N]    <- cft[N] / 2

      options$DIST$xMin    <- xMin
      options$DIST$xMax    <- xMax
      options$DIST$xMean   <- xMean
      options$DIST$cft     <- cft
      options$DIST$cfLimit <- cfLimit
    }


    # ALGORITHM ---------------------------------------------------------------

    if (length(x) == 0) {
      #if (missing(x)) {
      x <- seq(xMin, xMax, length.out = options$xN)
    }

    if (options$isInterp) {
      xOrg <- x
      #Chebysev points
      x <-
        (xMax - xMin) * (-cos(pi * (0:options$chebyPts) / options$chebyPts) + 1) / 2 + xMin
    } else {
      xOrg <- c()
    }

    #WARNING: OUT of range
    if (any(x < xMin) || any(x > xMax)) {
      warning(
        "CharFun: cf2DistGP" ,
        "x out of range (the used support): [xMin, xMax] = [",
        xMin,
        ", ",
        xMax,
        "]!"
      )
    }

    # Evaluate the required functions
    szx <- dim(x)
    x <- c(x)
    E <- exp((-1i) * x %*% t(t))

    # CDF estimate computed by using the simple trapezoidal quadrature rule
    cdf <- (xMean - x) / 2 + Im(E %*% (cft / t))
    cdf <- 0.5 - (cdf %*% dt) / pi

    # Correct the CDF (if the computed result is out of (0,1))
    # This is useful for circular distributions over intervals of length 2*pi,
    # as e.g. the von Mises distribution
    corrCDF <- 0
    if (options$correctedCDF) {
      if (min(cdf) < 0) {
        corrCDF <- min(cdf)
        cdf <- cdf - corrCDF
      }
      if (max(cdf) > 1) {
        corrCDF <- max(cdf) - 1
        cdf <- cdf - corrCDF
      }
    }

    dim(cdf) <- szx

    # PDF estimate computed by using the simple trapezoidal quadrature rule
    pdf <- 0.5 + Re(E %*% cft)
    pdf <- (pdf %*% dt) / pi
    pdf[pdf < 0] <- 0
    dim(pdf) <- szx
    dim(x) <- szx

    # REMARK:
    # Note that, exp(-1i*x_i*0) = cos(x_i*0) + 1i*sin(x_i*0) = 1. Moreover,
    # cf(0) = 1 and lim_{t -> 0} cf(t)/t = E(X) - x. Hence, the leading term of
    # the trapezoidal rule for computing the CDF integral is CDFfun_1 = (xMean- x)/2,
    # and PDFfun_1 = 1/2 for the PDF integral, respectively.

    # Reset the transformed CF, PDF, and CDF to the original values
    if (options$isCompound) {
      cf2  <- cfOld
      cdf <- cfLimit + cdf * (1 - cfLimit)
      pdf <- pdf * (1 - cfLimit)
      pdf[x == 0] = Inf
      pdf[x == xMax] = NaN
    }

    # Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
    # PrecisionCrit should be small for t > T, smaller than tolerance
    # options$crit
    PrecisionCrit <- abs(cft[length(cft)] / t[length(t)])
    isPrecisionOK <- (PrecisionCrit <= options$crit)

    # QF evaluated by the Newton-Raphson iterative scheme
    if (length(prob) > 0) {
      isPlot <- options$isPlot
      options$isPlot <- FALSE
      isInterp <- options$isInterp
      options$isInterp <- FALSE
      szp <- dim(prob)
      prob <- c(prob)
      maxiter <- options$maxiter
      crit <- options$crit
      qf <- options$qf0
      criterion <- TRUE
      count <- 0
      res <- cf2DistGP(cf2, qf, options = options)
      cdfQ <- res$cdf
      pdfQ <- res$pdf


      while (criterion) {
        count <- count + 1
        correction <- (cdfQ - corrCDF - prob) / pdfQ
        qf <- pmax(xMin, pmin(xMax, qf - correction))

        res <- cf2DistGP(function(x)
          cf2(x), x = qf, options = options)
        cdfQ <- res$cdf
        pdfQ <- res$pdf

        criterion <- any(abs(correction) > crit * abs(qf)) &&
          max(abs(correction)) > crit &&
          count < maxiter
      }

      dim(qf) <- szp
      dim(prob) <- szp
      options$isPlot <- isPlot
      options$isInterp <- isInterp

    } else {
      qf <- c()
      count = c()
      correction = c()
      #prob = c()
    }

    if (options$isInterp) {
      id <- is.finite(pdf)
      PDF <-
        function(xNew)
          pmax(0, interpBarycentric(x[id], pdf[id], xNew)[[2]])

      id <- is.finite(cdf)
      CDF <-
        function(xNew)
          pmax(0, pmin(1, interpBarycentric(x[id]), cdf[id], xNew)[[2]])

      QF <- function(prob)
        interpBarycentric(cdf[id], x[id], prob)[[2]]

      RND <- function(n)
        QF(runif(n))
    }

    if (length(xOrg) > 0) {
      x <- xOrg
      cdf <- CDF(x)
      pdf <- PDF(x)
    }

    # Reset the correct value for compound PDF at 0
    if (options$isCompound) {
      pdf[x == 0] <- Inf
    }

    ## Result
    result <- list(
      "Description"        = 'CDF/PDF/QF from the characteristic function CF',
      "x"                  = x,
      "cdf"                = cdf,
      "pdf"                = pdf,
      "prob"               = prob,
      "qf"                 = qf,
      "cf"                 = cfOld,
      "isCompound"         = options$isCompound,
      "isCircular"         = options$isCircular,
      "isInterp"           = options$isInterp,
      "SixSigmaRule"       = options$SixSigmaRule,
      "N"                  = N,
      "dt"                 = dt,
      "T"                  = t[length(t)],
      "PrecisionCrit"      = PrecisionCrit,
      "myPrecisionCrit"    = options$crit,
      "isPrecisionOK"      = isPrecisionOK,
      "xMean"              = xMean,
      "xStd"               = xStd,
      "xMin"               = xMin,
      "xMax"               = xMax,
      "cfLimit"            = cfLimit,
      "cdfAdjust"          = correction,
      "nNewtonRaphsonLoops" = count,
      "options"            = options
    )

    if (options$isInterp) {
      result$PDF <- PDF
      result$CDF <- CDF
      result$QF <- QF
      result$RND <- RND
    }

    end_time <- Sys.time()
    result$tictoc <- end_time - start_time

    # PLOT the PDF / CDF

    if (length(x) == 1)
      options$isPlot = FALSE

    if (options$isPlot) {
      plot(
        x = x,
        y = pdf,
        main = "PDF Specified by the CF",
        xlab = "x",
        ylab = "pdf",
        type = "l",
        lwd = 2,
        col = "blue"
      )
      grid(lty = "solid")

      plot(
        x = x,
        y = cdf,
        main = "CDF Specified by the CF",
        xlab = "x",
        ylab = "cdf",
        type = "l",
        lwd = 2,
        col = "blue"
      )
      grid(lty = "solid")
    }

    return(result)

  }
