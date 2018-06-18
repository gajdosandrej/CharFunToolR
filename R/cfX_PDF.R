#' @title
#' Computes the characteristic function of the continuos distribution defined by its PDF function.
#'
#' @description
#' \code{cfX_PDF(t, pdfFun, method, tol)} computes the characteristic function of the continuos
#'  distribution defined by its PDF function, \code{pdfFun = function(x) pdf(x)}, here we assume \code{x >= 0},
#'  computed for real vector argument \code{t}.
#'
#' @param t real vector, where the characteristic function CF(t) will be evaluated.
#' @param pdfFun function handle used as the PDF function with the argument \code{x}, for \code{x >= 0}.
#' However, \eqn{pdfFun(z)} should accept as an input value any complex matrix \eqn{z},
#' @param method set the method used to compute the characteristic function,
#' either \code{method = 'def'} (or the standard definition of CF by using the pdf) or \code{method = 'fit'}
#' (by using the half-space Fourier integral transform of the pdf).
#' Default value is \code{method = 'def'}.
#' @param tol relative tolerance parameter used in the built-in R numerical integration algorithm \code{integrate}.
#' Default value is \code{tol = 1e-6}.
#'
#' @details
#' \code{cfX_PDF} is based on the standard integral representation of the characteristic function
#' of the continuous distribution defined by its PDF (here PDF is represented by the function handle
#' \eqn{pdfFun(x)} defined for \eqn{x >= 0}). Then,
#' \deqn{CF(t) = Integral_0^inf exp(i*t*x) * pdfFun(x) dx}.
#' Alternatively, by using the half-space Fourier Integral Transformation (FIT),
#' which could improve the highly oscillatory behaviour of the integrand function, we get
#' \deqn{CF(t) = Integral_0^inf (i/t) * exp(-x) * pdfFun(i*x/t) dx.}
#' If we define the integrand as \eqn{funCF(t, x) = (i / t) * exp(-x) * pdfFun(i * x / t)},
#' then by using a stabilizing transformation from \eqn{[0,inf]} to \eqn{[0,1]},
#' we can evaluate the CF by the following (well behaved) integral:
#' \deqn{CF(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.}
#'
#' Selection of the proper method (standard definition or the integral transformation FIT)
#' depends on the distribution and the particular form of its PDF.
#' \code{cfX_PDF} evaluates this integral by using the R built in function \code{integrate},
#' with precission specified by tolerance tol (default value is \code{tol = 1e-6}).
#'
#' @return
#' Returns values of characteristic function of the continuos distribution defined by its PDF function .
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
#' [3] WITKOVSKY, V., WIMMER, G., DUBY, T. Computing the aggregate loss distribution
#' based on numerical inversion of the compound empirical characteristic function
#' of frequency and severity. ArXiv preprint, 2017, arXiv:1701.08299.
#'
#' [4] DUBY, T., WIMMER, G., WITKOVSKY, V.(2016). MATLAB toolbox CRM for computing distributions
#' of collective risk models. Working Paper. Journal of Statistical Software.
#' arXiv preprint arXiv:1801.02248, 2018.
#'
#' @family Continuous Probability Distribution
#'
#' @example R/Examples/example_cfX_PDF.R
#'
#' @export
#'
cfX_PDF <- function(t, pdfFun, method, tol) {

        ## CHECK THE INPUT PARAMETERS

        if(missing(t)) {
                stop("VW:cfX_PDF:TooFewInputs")
        }
        if(missing(pdfFun)) {
                pdfFun <- function(x) exp(-x)
        }
        if(missing(method)) {
                method <- "definition"
        }
        if(missing(tol)) {
                tol <- 1e-6
        }
        reltol <- tol

        ## EVALUATE THE CHARACTERISTIC FUNCTION: cfX_PDF(t,pdfFun)

        sz <- dim(t)
        t  <- c(t)
        cf <- numeric()
        if(length(sz) > 0) {
                cf <- matrix(1, sz[1], sz[2])
        } else {
                cf <- rep(1, length(t))
        }
        id <- t!=0

        if(method == "def" || method == "definition" || method == "standard" || method == "direct") {
                for(i in 1:length(t)) {
                        if(id[i]) {
                                re <- integrate(function(x) {as.numeric(Re(funCFdef(pdfFun, t[i], (x / (1 - x)) ^ 2) * (2 * x / (1 - x) ^ 3)))}, 0, 1, rel.tol = reltol, subdivisions = 2000)$value
                                im <- integrate(function(x) {as.numeric(Im(funCFdef(pdfFun, t[i], (x / (1 - x)) ^ 2) * (2 * x / (1 - x) ^ 3)))}, 0, 1, rel.tol = reltol, subdivisions = 2000)$value
                                cf[i] <- re + 1i*im
                        }
                }
        } else if(method == "fit" || method == "fourier" || method == "transoform") {
                for(i in 1:length(t)) {
                        if(id[i]) {
                                re <- integrate(function(x) {as.numeric(Re(funCFfit(pdfFun, t[i], (x / (1 - x)) ^ 2) * (2 * x / (1 - x) ^ 3)))}, 0, 1, rel.tol = reltol, subdivisions = 2000)$value
                                im <- integrate(function(x) {as.numeric(Im(funCFfit(pdfFun, t[i], (x / (1 - x)) ^ 2) * (2 * x / (1 - x) ^ 3)))}, 0, 1, rel.tol = reltol, subdivisions = 2000)$value
                                cf[i] <- re + 1i*im
                        }
                }
        } else {
                for(i in 1:length(t)) {
                        if(id[i]) {
                                re <- integrate(function(x) {as.numeric(Re(funCFdef(pdfFun, t[i], (x / (1 - x)) ^ 2) * (2 * x / (1 - x) ^ 3)))}, 0, 1, rel.tol = reltol, subdivisions = 2000)$value
                                im <- integrate(function(x) {as.numeric(Im(funCFdef(pdfFun, t[i], (x / (1 - x)) ^ 2) * (2 * x / (1 - x) ^ 3)))}, 0, 1, rel.tol = reltol, subdivisions = 2000)$value
                                cf[i] <- re + 1i*im                     }
                }
        }

        dim(cf) <- sz

        return(as.complex(cf))

}

## Function funCFdef
funCFdef <- function(pdfFun, t, x) {

        # funCFdef Integrand function of the integral representation of the
        # characteristic function CF defined by using the pdfFun(x),
        # for x >=0, and the real (vector) argument t.

        ## ALGORITHM

        x <- c(Conj(x))
        f <- t(pdfFun(x) * t(exp(1i * t %*% t(x))))

        return(f)
}

## Function funCFfit
funCFfit <- function(pdfFun, t, x) {

        # funCFfit Integrand function of the integral representation of the characteristic function CF
        # of the distribution defined by using the half-space Fourier Integral Transform (FIT)
        # of its pdfFun(x), for x >=0, and the real (vector) argument t.

        ## ALGORITHM

        ti <- 1 / c(t)
        x <- c(Conj(x))

        ot <- rep(1, length(ti))
        ox <- rep(1, length(x))
        f <- (1i * ti %*% t(ox)) * pdfFun(1i * ti %*% t(x)) * exp(-ot %*% t(x))

        return(f)
}
