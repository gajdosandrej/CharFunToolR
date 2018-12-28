#' @title Characteristic function of the weighted mixture distribution
#' of independent DIRAC random variables concentrated at the fixed constants
#'
#' @description
#' Characteristic function of the weighted mixture distribution
#' of independent DIRAC random variables \eqn{D_1,...,D_N}, concentrated
#' at the fixed constants (data) given by the vector \eqn{d = (d_1,...,d_N)}.
#'
#' That is, \eqn{cf(t) = weight_1*cfD(d_1*t) +...+ weight_N*cfD(d_N*t)},
#' where \eqn{cfD(t)} represents the characteristic function of the DIRAC RV
#' concentrated at the constant \eqn{d=1}, i.e. \eqn{cfD(t) = exp(1i*t)}.
#'
#' \code{cfE_DiracMixture(t, d, weight, cfX)} evaluates the compound characteristic function
#' \deqn{cf(t) = cfE_DiracMixture(-1i*log(cfX(t)),d,weight) = weight_1*cfX(t)^d_1 +...+ weight_N*cfX(t)^d_N}
#' where \code{cfX} denotes the function handle of the characteristic function \eqn{cfX(t)} of the random variable \eqn{X}.
#'
#' @family Empirical Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Empirical_characteristic_function}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param d vector of constants (data) where the DIRAC RVs are concentrated.
#' If empty, default value is \code{d = 1}.
#' @param weight vector of weights of the distribution mixture.
#' If empty, default value is \code{weight = 1/length(d)}.
#' @param cfX function handle of the characteristic function of a random
#' variable \eqn{X}. If \code{cfX} is non-empty, a compound CF is evaluated as \eqn{cf(t) = cf(-1i*log(cfX(t)),d,weight)}.
#'
#' @return Characteristic function \eqn{cf(t)} of the EMPIRICAL distribution, based on the observed data.
#'
#' @references
#' WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
#' loss distribution based on numerical inversion of the compound empirical
#' characteristic function of frequency and severity. arXiv preprint arXiv:1701.08299.
#'
#' @note Ver.: 23-Sep-2018 14:59:46 (consistent with Matlab CharFunTool v1.3.0, 23-Jun-2017 10:00:49).
#'
#' @example R/Examples/example_cfE_DiracMixture.R
#'
#' @export
#'
cfE_DiracMixture <- function(t, d, weight, cfX) {
        ## CHECK THE INPUT PARAMETERS
        if (missing(d)) {
                d <- vector()
        }
        if (missing(weight)) {
                weight <- vector()
        }
        if(missing(cfX)) {
                cfX <- NULL
        }

        ##
        if (length(d) == 0) {
          d <- 1
        }
        if(length(weight) == 0) {
                weight <- 1 / length(d)
        }
        if(length(weight) == 1) {
                d <- Conj(d)
        } else {
                l_max <- max(c(length(d), length(weight)))
                if (l_max > 1) {
                        if (length(d) == 1) {
                                d <- rep(d, l_max)
                        }
                        if (length(weight) == 1) {
                                weight <- rep(weight, l_max)
                        }
                        if (any(lengths(list(d, weight)) < l_max)) {
                                stop("Input size mismatch.")
                        }
                }
        }

        # Special treatment for mixtures with large number of variables

        szcoefs <- dim(t(d))
        szcoefs <- szcoefs[1] * szcoefs[2]

        szt <- dim(t)
        sz <- dim(t(t))[1] * dim(t(t))[2]

        szcLimit <- ceiling(1e3 / (sz / 2 ^ 16))
        idc <- (1:(trunc(szcoefs / szcLimit) + 1))

        # Characteristic function of a weighted mixture of Dirac variables

        t <- c(t)
        idx0 <- 1
        cf <- 0

        for (j in idc) {
                idx1 <- min(idc[j] * szcLimit, szcoefs)
                idx <- (idx0:idx1)
                idx0 <- idx1 + 1

                if (is.null(cfX)) {
                        aux <- exp(1i * t %*% t(d[idx]))
                } else {
                        aux <- t(apply(as.matrix(cfX(t)), 1, FUN = '^', d[idx]))
                }

                if (length(weight) == 1) {
                        cf <- cf + apply(weight * aux, 1, sum)
                } else {
                        cf <- cf + apply(aux * weight[idx], 1, sum)
                        cf <- cf + apply(t(apply(aux * weight[idx], 1, FUN = '*', c(1,2,3))), 1, sum)
                }

        }

        dim(cf) <- szt

        return(cf)
}
