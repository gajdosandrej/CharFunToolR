#' @title
#' Evaluating CDF/PDF/QF (quantiles) from the characteristic function CF
#' of a (continuous) distribution F by using the Fast Fourier Transform (FFT) algorithm
#'
#' @description
#' TEST VERSION !
#' \code{cf2DistFFT(cf, x, prob, options)} evaluates the approximate values CDF(x), PDF(x),
#' and/or the quantiles QF(prob) for given \code{x} and \code{prob}, by interpolation
#' from the PDF-estimate computed by the numerical inversion of the given
#' characteristic function CF by using the FFT algorithm.
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
#' @param x       vector of values from domain of the distribution \eqn{F}, if \code{x} is left empty,
#' \code{cf2DistFFT} automatically selects vector \code{x} from the domain.
#' @param prob    vector of values from \eqn{[0,1]} for which the quantiles
#' will be estimated, if \code{prob} is left empty, \code{cf2DistFFT} automatically
#' selects vector \code{prob = c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999)}.
#' @param options  structure (list) with the following default parameters:
#' \itemize{
#'     \item \code{options$isCompound = FALSE} treat the compound distributions, of the RV \eqn{Y = X_1 + ... + X_N},
#'     where \eqn{N} is discrete RV and \eqn{X\ge 0} are iid RVs from nonnegative continuous distribution,
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
#'         \item \code{options$DIST$xMin = xMin} the lower limit of \code{x},
#'         \item \code{options$DIST$xMax = XMax} the upper limit of \code{x},
#'         \item \code{options$DIST$xMean = xMean} the MEAN value of \code{x},
#'         \item \code{options$DIST$cft = cft} CF evaluated at \eqn{t_j} : \eqn{cf(t_j)}.
#'         }
#'     }
#'
#' @return
#' \item{result}{structure (list) with with the following results values:}
#' \item{result$x = x;}{}
#' \item{result$cdf = cdf;}{}
#' \item{result$pdf = pdf;}{}
#' \item{result$prob = prob;}{}
#' \item{result$qf = qf;}{}
#' \item{result$xFTT = xFFT;}{}
#' \item{result$pdfFFT = pdfFFT;}{}
#' \item{result$cdfFFT = cdfFFT;}{}
#' \item{result$SixSigmaRule = options.SixSigmaRule;}{}
#' \item{result$N = N;}{}
#' \item{result$dt = dt;}{}
#' \item{result$T = t[length(t)];}{}
#' \item{result$PrecisionCrit = PrecisionCrit;}{}
#' \item{result$myPrecisionCrit = options.crit;}{}
#' \item{result$isPrecisionOK = isPrecisionOK;}{}
#' \item{result$xMean = xMean;}{}
#' \item{result$xStd = xStd;}{}
#' \item{result$xMin = xMin;}{}
#' \item{result$xMax = xMax;}{}
#' \item{result$cf = cf;}{}
#' \item{result$options = options;}{}
#' \item{result$tictoc = toc.}{}
#'
#' @details
#' The outputs of the algorithm \code{cf2DistFFT} are approximate values!
#' The precission of the presented results depends on several different factors: \cr
#' - application of the FFT algorithm for numerical inversion of the CF \cr
#' - selected number of points used by the FFT algorithm (by default \code{options$N = 2^10}), \cr
#' - estimated/calculated domain \eqn{[A,B]} of the distribution \eqn{F}, \cr
#' used with the FFT algorithm. Optimally, \eqn{[A,B]} covers large part of the distribution domain,
#' say more than \eqn{99\%}. However, the default automatic procedure for selection of the domain
#' \eqn{[A,B]} may fail. It is based on the 'SixSigmaRule': \eqn{A = MEAN - SixSigmaRule * STD},
#' and \eqn{B = MEAN + SixSigmaRule * STD}. Alternatively, change the \code{options$SixSigmaRule}
#' to different value, say \eqn{12}, or use the \code{options$xMin} and \code{options$xMax} to set
#' manually the values of \eqn{A} and \eqn{B}.
#'
#' @references
#' [1] WITKOVSKY, V.: On the exact computation of the density and
#' of the quantiles of linear combinations of t and F random variables.
#' Journal of Statistical Planning and Inference 94 (2001), 1-13.
#'
#' [2] WITKOVSKY, V.: Matlab algorithm TDIST: The distribution of a linear combination
#' of Student's t random variables. In COMPSTAT 2004 Symposium (2004), J. Antoch, Ed.,
#' Physica-Verlag/Springer 2004, Heidelberg, Germany, pp. 1995-2002.
#'
#' [3] WITKOVSKY, V., WIMMER, G., DUBY, T.: Logarithmic Lambert W x F random variables
#' for the family of chi-squared distributions and their applications.
#' Statistics & Probability Letters 96 (2015), 223-231.
#'
#' [4] WITKOVSKY, V. (2016): Numerical inversion of a chracteristic function:
#' An alternative tool to form the probability distribution of output quantity
#' in linear measurement models. Acta IMEKO, 5(3), 32-34.
#'
#' [5] WITKOVSKY, V., WIMMER G., DUBY T. (2016): Computing the aggregate loss distribution
#' based on numerical inversion of the compound empirical characteristic function
#' of frequency and severity. Preprint submitted to Insurance: Mathematics and Economics.
#'
#' [6] DUBY, T., WIMMER, G., WITKOVSKY, V. (2016): MATLAB toolbox CRM for computing distributions
#' of collective risk models. Preprint submitted to Journal of Statistical Software.
#'
#' @example R/Examples/example_cf2DistFFT.R
#'
#' @export
#'
cf2DistFFT <- function(cf, x, prob, options) {
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
                options$tolDiff = 1e-4
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
        if (is.null(options$isPlotFFT)) {
                options$isPlotFFT <- FALSE
        }

        if (is.null(options$xN)) {
                options$xN <- 101
        }


        if (is.null(options$chebyPts)) {
                options$chebyPts <- 2^9
        }


        if (is.null(options$isInterp)) {
                options$isInterp <- FALSE
        }


        ## GET/SET the DEFAULT parameters and the OPTIONS
        # First, set a special treatment if the real value of CF at infinity (largevalue)
        # is positive, i.e. const = real(cf(Inf)) > 0. In this case the compound CDF
        # has jump at 0 of size equal to this value, i.e. cdf(0) = const, and pdf(0) = Inf.
        # In order to simplify the calculations, here we calculate PDF and CDF of a distribution
        # given by transformed CF, i.e. cf_new(t) = (cf(t)-const) / (1-const);
        # which is converging to 0 at Inf, i.e. cf_new(Inf) = 0. Using the transformed CF
        # requires subsequent recalculation of the computed CDF and PDF, in order to get
        # the originaly required values: Set pdf_original(0) =  Inf &
        # pdf_original(x) = pdf_new(x) * (1-const), for x > 0.
        # Set cdf_original(x) =  const + cdf_new(x) * (1-const).
        #
        const <- Re(cf(1e30))
        cfOld <- NULL

        if (options$isCompound) {
                cfOld <- cf
                if (const > 1e-13) {
                        cf2 <- cf
                        cf <- function(t) {
                                (cf2(t) - const) / (1 - const)
                        }
                }
        }

        if (length(options$DIST) > 0) {
                xMean <- options$DIST$xMean
                cft <- options$DIST$cft
                xMin <- options$DIST$xMin
                xMax <- options$DIST$xMax
                N <- length(cft)
                k <- 0:(N-1)
                xRange <- xMax - xMin
                dt <- 2 * pi / xRange
                t <- (k - N/2 + 0.5) * dt
                xStd <- numeric()
        } else {
                N <- 2*options$N
                dt <- options$dt
                T <- options$T
                xMin <- options$xMin
                xMax <- options$xMax
                xMean <- options$xMean
                xStd <- options$xStd
                SixSigmaRule <- options$SixSigmaRule
                tolDiff <- options$tolDiff
                cft <- cf(tolDiff*(1:4))

                if (length(xMean) == 0) {
                        xMean <- Re((-cft[2] + 8*cft[1]-8*Conj(cft[1]) + Conj(cft[2]))/(1i*12*tolDiff))
                }
                if (length(xStd) == 0) {
                        xM2 <- Re(-(Conj(cft[4]) - 16*Conj(cft[3]) + 64*Conj(cft[2])
                                      + 16*Conj(cft[1]) - 130 + 16*cft[1] + 64*cft[2]
                                      - 16*cft[3]+cft[4])/(144*tolDiff^2))
                        xStd  <- sqrt(xM2 - xMean^2)
                }
                if (is.finite(xMin) && is.finite(xMax)) {
                        xRange <- xMax - xMin
                } else if (length(T) > 0) {
                        xRange <- 2 * pi / (2 * T / N)
                        if (is.finite(xMax)) {
                                xMin <- xMax - xRange
                        } else if (isfinite(xMin)) {
                                xMax <- xMin + xRange
                        } else {
                                xMin <- xMean - xRange/2
                                xMax <- xMean + xRange/2
                        }
                } else if (length(dt) > 0) {
                        xRange <- 2 * pi / dt
                        if (is.finite(xMax)) {
                                xMin <- xMax - xRange
                        } else if (isfinite(xMin)) {
                                xMax <- xMin + xRange
                        } else {
                                xMin <- xMean - xRange/2
                                xMax <- xMean + xRange/2
                        }
                } else {
                        if (is.finite(xMin)) {
                                xMax <- xMean + SixSigmaRule * xStd
                        } else if (is.finite(xMax)) {
                                xMin <- xMean - SixSigmaRule * xStd
                        } else {
                                xMin <- xMean - SixSigmaRule * xStd
                                xMax <- xMean + SixSigmaRule * xStd
                        }
                        xRange <- xMax - xMin
                }

                dt      <- 2 * pi / xRange
                k       <- 0:(N-1)
                t       <- (k - N/2 + 0.5) * dt
                cft     <- cf(t[(N/2+1):length(t)])
                cft     <- c(Conj(cft[seq(length(cft), 1, -1)]), cft)
                #cft[1]  <- cft[1] / 2
                #cft[N]  <- cft[N] / 2

                options$DIST$xMin    <- xMin
                options$DIST$xMax    <- xMax
                options$DIST$xMean   <- xMean
                options$DIST$cft     <- cft
        }


        # ALGORITHM ---------------------------------------------------------------

        A <- xMin
        B <- xMax

        dx <- (B - A) / N
        c <- (as.complex(-1))^(A * (N - 1) / (B - A)) / (B-A)

        C <- c * (as.complex(-1))^((1 - 1 / N) * k)
        D <- (as.complex(-1))^(-2 * (A / (B - A)) * k)
        pdfFFT <- pmax(0, Re(C * stats::fft(D * cft)))
        cdfFFT <- pmin(1, pmax(0, 0.5 + Re(1i * C * stats::fft(D * cft / t))))


        xFFT <- A + k * dx

        # Reset the transformed CF, PDF, and CDF to the original values
        if(options$isCompound) {
                cf <- cfOld
                cdfFFT <- const + cdfFFT * (1 - const)
                pdfFFT <- pdfFFT * (1 - const)
                pdfFFT[x==0] <- Inf
        }
      #  print(cdfFFT)

        # Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
        # PrecisionCrit should be small for t > T, smaller than tolerance
        # options.crit
        PrecisionCrit <- abs(cft[length(cft)] / t[length(t)])
        isPrecisionOK <- (PrecisionCrit<=options$crit)

        ## INTERPOLATE QUANTILE FUNCTION for required prob values: QF(prob)
        if(length(prob)==0) {
                prob <- c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999)
        }



        cdfU <- sort(unique(cdfFFT))

       # print(cdfU)

        id <- firstOccIdx(cdfFFT)
        xxU <- xFFT[id]
        szp <- dim(prob)
        eps <- 2.2204e-16



        qfFun <- function(prob) {pracma::pchip(c(-eps, cdfU), c(-eps, xxU+dx/2), prob)}
        qf <- qfFun(prob)
        dim(qf) <- szp

        # In the next two interpolations, pchip interpolation method (function) is used
        # instead of (default) linear interpolation

        # INTERPOLATE CDF required values: CDF(x)
        if(length(x)==0) {
                xempty<-TRUE
                x <- seq(from = xMin, to = xMax, length.out = options$xN)
        }
        else{
                xempty<-FALSE
        }


        if (is.null(options$isInterp)) {
                x0<-x
                # Chebyshev points
                x<-(xMax - xMin)*(-cos(pi*seq(0,options$chebyPts)/options$chebyPts)+1)/2+xMin

        }
        else{
            x0<-vector()
        }


        # WARNING: OUT - of - range
        if(any(xMin)||any(x>xMax)){
                error = function(err) {print(paste("cf2DistFFT: x out-of-range:
        [xMin, xMax] = [",xMin,", ",xMax,"] !", err))
                }
        }



        szx <- dim(x)
        # cdfFun = @(x) interp1([-eps;xxU+dx/2],[-eps;cdfU],x(:));
        cdfFun <- function(x) {pracma::pchip(sort(xxU), cdfU, c(x))}
        cdf <- cdfFun(x)
        dim(cdf) <- szx

        # TRY INTERPOLATE PDF required values: PDF(x)
        tryCatch({pdfFun <- function(x) {pracma::pchip(sort(xFFT), pdfFFT, c(x))}
                 pdf <- pmax(0, pdfFun(x))
                 dim(pdf) <- szx},
                 error = function(err) {print(paste("cf2DistFFT: Unable to interpolate the required PDF values", err))
                         pdf <- NaN*x})
        dim(x) <- szx



        if (options$isInterp) {
                id<-is.finite(pdf)
                PDF<-function(xnew){pracma::pchip(sort(xnew),x[id],pdf[id])}
                id<-is.finite(cdf)
                CDF<-function(xnew){pracma::pchip(sort(xnew),x[id],cdf[id])}
                QF<-function(prob){pracma::pchip(sort(prob),x[id],cdf[id])}
                RND<-function(dim){pracma::pchip(sort(dim),x[id],cdf[id])}

                tryCatch({if(length(x)>0){
                 x<-x0
                 cdf<-CDF[x]
                 pdf<-PDF[x]
                }
                },
                error=function(err){print(paste("cf2DistGPT: Problem using the interpolant function",err))})
        }
        else{
         PDF<-vector()
         CDF<-vector()
         QF<-vector()
         RND<-vector()
        }

        ## Result
        result <- list(
                "Description"="CDF/PDF/QF from the characteristic function CF",
                "inversionMethod"="Discrete version of standard inversion formula",
                "quadratureMethod"="Fast Fourier Transform (FFT) algorithm",
                "x"                  = x,
                "cdf"                = cdf,
                "pdf"                = pdf,
                "prob"               = prob,
                "qf"                 = qf,

                if(is.null(options$isInterp)){
                        "PDF" =  PDF
                        "CDF" = CDF
                        "QF"=QF
                        "RND"=RND
                },
                "xFFT"               = xFFT,
                "pdfFFT"             = pdfFFT,
                "cdfFFT"             = cdfFFT,
                "SixSigmaRule"       = options$SixSigmaRule,
                "N"                  = options$N,
                "dt"                 = dt,
                "T"                  = t[length(t)],
                #"PrecisionCrit"      = PrecisionCrit,
                #"myPrecisionCrit"    = options$crit,
                "isPrecisionOK"      = isPrecisionOK,
                "xMean"              = xMean,
                "xStd"               = xStd,
                "xMin"               = xMin,
                "xMax"               = xMax,
                "cf"                 = cf,
                "const"              = const,
                "isCompound"         = options$isCompound,
                "options"            = options

                  )

        end_time <- Sys.time()
        result$tictoc <- end_time - start_time

        ## PLOT THE PDF/CDF, if required

        if (length(x) == 1)
                options$isPlot = FALSE

        if(options$isPlotFFT) {
                x <- xFFT
                pdf <- pdfFFT
                cdf <- cdfFFT
        }

        if (options$isPlot) {
                plot(
                        x = x,
                        y = pdf,
                        main = "PDF Specified by the Characteristic Function CF",
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
                        main = "CDF Specified by the Characteristic Function CF",
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



firstOccIdx <- function(x) {
        x_uniqe <- sort(unique(x))
        indices <- vector()
        for (i in 1:length(x_uniqe)) {
                idx <- 1
                for (j in 1:length(x)) {
                        if (x[j] == x_uniqe[i]) {
                                idx <- j
                                indices <- c(indices, idx)
                                break
                        }
                }
        }
        return(indices)
}

