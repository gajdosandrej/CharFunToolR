#' @title Characteristic function of the OGIVE EMPIRICAL distribution based on the observed histogram
#'
#' @description
#' Characteristic function of the OGIVE EMPIRICAL distribution based
#' on the observed histogram (given by bins and frequencies).
#' That is, \eqn{cf(t)} is given as a weighted mixture of the UNIFORM iid RVs.
#' \deqn{cf(t) = cfE_EmpiricalOgive(t,bins,freq) = sum_{j=1}^n freq_j * cf_Uniform(t,bins_{j-1},bins_{j}),}
#' where \eqn{cf_Uniform(t,bin_{j-1},bin_{j})} is CF of the Uniform distribution on the interval \eqn{[bins_{j-1},bins_j]}.
#'
#' The bins and frequencies = counts/sum(counts) are based on the discretized/grouped
#' or histogram data, for more details see e.g. the R function \code{hist}.
#'
#' @family Empirical Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Empirical_characteristic_function}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param bins vector of breaks of histogram bins.
#' @param freq frequencies for each bin of histogram.
#'
#' @return Characteristic function \eqn{cf(t)} of the EMPIRICAL distribution, based on the observed data.
#'
#' @references
#' WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
#' loss distribution based on numerical inversion of the compound empirical
#' characteristic function of frequency and severity. arXiv preprint arXiv:1701.08299.
#'
#' @note Ver.: 02-Oct-2018 14:07:43 (consistent with Matlab CharFunTool v1.3.0, 24-Jun-2017 09:34:15).
#'
#' @example R/Examples/example_cfE_EmpiricalOgive.R
#'
#' @export
#'
cfE_EmpiricalOgive <- function(t, bins, freq) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(bins)) {
                stop("bins argument must be nonempty!")
        }
        if(missing(freq)) {
                stop("freq argument must be nonempty!")
        }
        if((length(bins) - 1) != length(freq)) {
                stop("Input size mismatch.")
        }

        ## Characteristic function
        szt <- dim(t)
        t <- c(t)
        cf <- 0

        n <- length(bins)
        n1 <- n-1

        aux <- exp(1i * t %*% t(bins))
        aux <- (aux[,2:n] - aux[,1:n1]) / (1i * t %*% t((bins[2:n] - bins[1:n1])))
        if(length(freq) == 1) {
                cf <- cf + colSums(freq * aux)
        } else {
                cf   = cf + apply(t(t(aux) * freq), 1, "sum")
        }

        cf[t==0] <- 1
        dim(cf) <- szt

        return(cf)
}
