#' @title
#' Characteristic function of a linear combination (convolution) of iid random variables X
#'
#' @description
#' \code{cf_Conv(t, cfX, coef, n)} computes the characteristic function of a linear combination
#' (convolution) of iid random variables \eqn{X} with given (common) characteristic function \eqn{cfX(t)}
#' and the coefficients coef, and such that \eqn{Y = coef(1)*X + ... + coef(N)*X},
#' i.e. \eqn{cf = cfX(coef(1)*t) * ... * cfX(coef(N)*t)}.
#'
#' If or all \code{coef == 1}, then \eqn{cf = cfX(t)^N}, . Moreover, with using the optional parameter \code{n},
#' the algorithm evaluates CF of the convolution of \eqn{n} independent copies of the random variable \eqn{Y},
#' i.e. \eqn{Z = Y + ... + Y}, where \eqn{Y = coef(1)*X + ... + coef(N)*X},
#' i.e. \eqn{cf = (cfX(coef(1)*t) * ... * cfX(coef(N)*t))^n}.
#'
#' @param t vector (or array) of input values \code{t} where the \code{cf_Conv} is evaluated.
#' @param cfX function handle to the given chracteristic function \eqn{cfX(t)}.
#' @param coef vector of coeficients of the linear combination of the iid random variables,
#' such that \eqn{Y = sum(coef(k) * X_k), cf = Prod_{k=1}^N cfX(coef(k)*t)}.
#' @param n optional power coeficient of additional convolution of the combiend CF of \eqn{Y}.
#' With using \eqn{n} (if empty, default value is \code{n = 1}), \eqn{cf = (Prod_{k=1}^N cfX(coef(k)*t)^n}.
#'
#' @return
#' The characteristic function of a linear combination of iid RVs with characteristic function \eqn{cf},
#' evaluated at \eqn{t}.
#'
#' @family CF Tool
#'
#' @note Ver.: 06-Oct-2018 18:14:22 (consistent with Matlab CharFunTool v1.3.0, 9-May-2017 10:22:48).
#'
#' @example R/Examples/example_cf_Conv.R
#'
#' @export
#'
cf_Conv <- function(t, cfX, coef, n) {

        ## CHECK THE INPUT PARAMETERS
        if(missing(coef)) {
                coef <- vector()
        }
        if(missing(n)) {
                n <- numeric()
        }

        if(length(coef) == 0) {
                coef <- 1
        }


        ## Find the unique coefficients and their multiplicities
        nums <- numeric()
        if(length(coef) > 0 ) {
                coef <- sort(coef)
                coef_orig_sort <- coef
                m <- length(coef)
                coef <- unique(coef)
                idx <- firstOccIdx(coef_orig_sort)
                nums <- diff(c(idx,m+1))
        } else {
                nums <- 1
        }

        l_max <- max(c(length(coef), length(nums)))
        if (l_max > 1) {
                if (length(nums) == 1) {
                        nums <- rep(nums, l_max)
                }
                if (length(coef) == 1) {
                        coef <- rep(coef, l_max)
                }
                if ((any(lengths(list(coef, nums)) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        # Special treatment for linear combinations with large number of RVs
        szcoefs <- dim(as.matrix(coef))
        szcoefs <- szcoefs[1] * szcoefs[2]
        szt <- dim(as.matrix(t))
        sz <- szt[1] * szt[2]
        szcLimit <- ceiling(1e3 / (sz / (2 ^ 16)))
        idc <- 1:(as.integer(szcoefs / szcLimit) + 1)

        ## ALGORITHM
        szt <- dim(t)
        t <- c(t)
        idx0 <- 1
        cf <- 1
        for (j in 1:idc[length(idc)]) {
                idx1 <- min(idc[j] * szcLimit, szcoefs)
                idx <- idx0:idx1
                idx0 <- idx1 + 1
                aux <- t %*% t(coef[idx])
                aux <- t(t(cfX(aux)) ^ nums[idx])
                cf <- cf * apply(aux, 1, prod)
        }
        dim(cf) <- szt
        cf[t == 0] <- 1

        if (length(n) > 0) {
                if (length(n) == 1) {
                        cf <- cf ^ n
                } else {
                        stop("n should be a scalar (positive integer) value.")
                }
        }

        return(cf)
}

# auxiliary function to find the index of the first occurence of certain element (sorted sequence)
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



