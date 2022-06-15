#' @title Evaluates the PMF and CDF of a discrete random variable with lattice distribution
#'
#' @description
#' The \code{cf2PMF_FFT(cf, xMin, xMAX, xDelta, options)}  evaluates the distribution functions (PMF and CDF) of a
#' DISCRETE random variable (DRV) with lattice distribution defined on the
#' Support = \code{seq(xMin,xMax,xDelta)}, and fully specified by its characteristic
#' function (CF). PMF and CDF is evaluated at the specified support points
#' x = \code{seq(xMin,xMax,xDelta)} by numerical inversion of CF, based on using the
#' inverse FFT algorithm.. Algorithm considers only finite discrete support, so it
#' should be specified carefully and correctly (i.e., it should include a
#' probability mass equal to 1 - epsilon, where epsilon is an accepted
#' truncation error). See Warr (2014) for more information on the method used.
#'
#' @family CF Inversion Algorithm
#'
#' @importFrom stats runif
#' @importFrom graphics plot grid
#'
#' @seealso For more details see:
#' \url{https://arxiv.org/pdf/1701.08299.pdf}
#'
#' @param cf function handle of the characteristic function of the
#' discrete lattice distribution, defined on a subset of
#' integers.
#' @param xMin minimum value of the support of the (lattice) discrete RV X.
#' @param xMax maximum value (finite number) of the support of the (lattice)
#' DRV X.
#' @param xDelta \eqn{xDelta > 0} is the minimum difference of the support values
#' of the discrete RV X with lattice distribution given as Supp = \code{seq(xMin,xMax,xDelta)}.
#' If empty, the default value is \eqn{xDelta = 1}.
#' @param options structure with the following default parameters:
#' \itemize{
#' \item \code{options$xMin = 0}, minimum value of the support of X,
#' \item \code{options$xMax = 100}, maximum value of the support of X,
#' \item\code{options$xDelta = 1}, minimum difference of the support values of X,
#' \item \code{options$isPlot = true},   logical indicator for plotting PMF and CDF.
#' }
#'
#' @return
#' \itemize{
#' \item \code{result}  structure with all details,
#' \item \code{pmf}   vector of the PMF values evaluated at \code{x=seq(xMin, xMax)},
#' \item \code{cdf}   vector of the PMF values evaluated at \code{x=seq(xMin, xMax)}.
#' }
#'
#' @references
#' [1] Warr, Richard L. Numerical approximation of probability mass functions
#' via the inverse discrete Fourier transform. Methodology and Computing in
#' Applied Probability 16, no. 4 2014: 1025-1038.
#'
#' @example R/Examples/example_cf2PMF_FFT.R
#'
#' @export
#'
cf2PMF_FFT <- function(cf, xMin, xMax, xDelta,options) {

  ## CHECK/SET THE INPUT PARAMETERS

  start_time <- Sys.time()

  if (missing(xDelta)) {
    xDelta <- vector()
  }

  if (missing(xMax)) {
    xMax <- vector()
  }

  if (missing(xMin)) {
    xMin <- vector()
  }


  if (missing(options)) {
    options <- list()
  }

  if (is.null(options$xMin)) {
    options$xMin <- 0
  }
  if (is.null(options$xMax)) {
    options$xMax <- 100
  }

  if (is.null(options$xDelta)) {
    options$xDelta <- 1
  }

  if (is.null(options$isPlot)) {
    options$isPlot <- TRUE
  }

  if(length(xMax)==0){
    xMax<-options$xMax

  }

  if(length(xMin)==0){
    xMin<-options$xMin

  }

  if(length(xDelta)==0){
  xDelta<-options$xDelta

  }
  ## ALGORITHM
  x<-t(seq(xMin,xMax,1))
    N<-length(x)

    if(xMin==0){
      omega<-(x/xDelta)/N
      DFTfun<-cf(-2*pi*omega/xDelta)
    }

    else{
      xx<-x-xMin
      omega<-(xx/xDelta)/N
      DFTfun<-cf(-2*pi*omega/xDelta)*exp(1i*2*pi*omega*xMin/xDelta)
    }

    # PMF
    DFTfun<-c(DFTfun)
    eps<-2.2204e-16
    pmf   <- Re(pracma::ifft(DFTfun))
    pmf[pmf<100*eps]<-0

    # CDF
    cdf    <- cumsum(pmf)
    tictoc<-Sys.time()-start_time

## RESULTS
    result <- list(
      "Description"="PMF/CDF of a discrete distribution from its CF",
      "inversionMethod"="inverse FFT algorithm",
      "pmf"="pmf",
      "cdf"="cdf",
      "x"="x",
      "cf"="cf",
      "DFTfun"="DFTfun",
      "options"="options",
      "tictoc"="tictoc")

      ## PLOT THE PDF/CDF
    if (length(x) == 1){
      options$isPlot = FALSE
    }
    if (options$isPlot) {
      # PMF
      plot(
        x = x,
        y = pmf,
        main = "PMF Specified by the Characteristic Function CF",
        xlab = "x",
        ylab = "pmf",
        type = "l",
        lwd = 2,
        col = "blue"

)

    # CDF
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
    }
    return(result)
}
