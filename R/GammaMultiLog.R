#' @title
#' Evaluates logarithm of the multivariate gamma function of order p
#'
#' @description
#' \code{GammaMultiLog(z, p, funmode)} evaluates logarithm of the multivariate gamma function
#' of order \code{p}, the argument \code{z} may be complex and of any size.
#'
#' Multivariate gamma function of order \eqn{p}, say \eqn{gamma_p(z)},
#' is defined as \eqn{gamma_p(z) = pi^(p*(p-1)/2) * prod_{j=1}^p gamma(z-(j-1)/2)}.
#' Hence, the logarithm of multivariate gamma
#' is given by \eqn{log(gamma_p(z))=(p*(p-1)/2)*log(pi) + sum_{j=1}^p log(gamma(z-(j-1)/2))}.
#'
#' @family Utility Function
#'
#' @param z complex argument, of any size (vector, matrix).
#' @param p order of the multivariate gamma, if empty, default value is \code{p = 1}.
#' @param funmode function mode: If \code{funmode = 0} \eqn{f = log[gamma_p(z))}. If \code{funmode = 1} \eqn{f = gamma_p(z)}.
#'
#' @return  Function returns values of the logarithm of multivariate gamma function of order \code{p} evaluated in points \code{z}.
#'
#' @note Ver.: 18-Sep-2019 17:14:16 (consistent with Matlab CharFunTool v1.3.0, 17-Aug-2018 19:45:37).
#'
#' @example R/Examples/example_GammaMultiLog.R
#'
#' @export
#'
GammaMultiLog <- function(z, p, funmode) {
        ## CHECK THE INPUT PARAMETERS

        if(missing(funmode)) {
                funmode <- numeric()
        }

        if(missing(p)) {
                p <- numeric()
        }

        if(length(p) == 0) {
                p <- 1
        }

        if(length(funmode) == 0) {
                funmode <- 0
        }

        sz <- dim(z)
        z <- c(z)

        # This was corrected on September 18, 2019 (from f = (p*(p-1)/2)*log(pi))
        f <- (p * (p - 1) / 4) * log(pi)

        for(j in 1:p) {
                f <- f + GammaLog(z - (j - 1) / 2)
        }

        if(funmode == 1) {
                f <- exp(f)
        }

        dim(f) <- sz

        return(f)
}
