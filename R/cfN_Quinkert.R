#' @title Characteristic function of the Quinkert distribution
#'
#' @description
#'\code{cfN_Quinkert(t,a,b)} evaluates the characteristic function \eqn{cf(t)}
#' of the Quinkert distribution, with the parameters \code{a} (a > 0) and \code{b} (b > 0), i.e.
#'
#'  \eqn{cf(t) = cfN_Quinkert(t,a,b) = 1F1(a,a+b,e^(1i*t)-1)}, where \eqn{1F1} denotes the confluent hypergeometric (Kummer's) function.
#'   For more details see [4], p. 564.
#'
#' @family Discrete Probability Distribution
#'
#' @references
#'  [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
#'   distribution based on numerical inversion of the compound empirical
#'    characteristic function of frequency and severity. Preprint submitted
#'    to Insurance: Mathematics and Economics.
#'
#'[2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
#'    computing distributions of collective risk models. Preprint submitted
#'    to Journal of Statistical Software.
#'
#'[3] WITKOVSKY V. (2016). Numerical inversion of a characteristic
#'    function: An alternative tool to form the probability distribution of
#'    output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#'[4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
#'    probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN 3-87773-025-6.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param a vector of the 'shape' parameters \code{a > 0}. If empty, default value is \code{a = 1}.
#' @param b vector of the 'shape' parameters \code{b > 0}. If empty, default value is \code{b = 1}.
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Quinkert distribution.
#'
#' @note ver.: 31-Jul-2021 12:36:48 (consistent with Matlab CharFunTool v1.5.1, 15-Nov-2016 13:36:26).
#'
#' @example R/Examples/example_cfN_Quinkert.R
#'
#' @export
#'
cfN_Quinkert <- function(t, a = 1, b = 1, cfX) {
        ## CHECK THE INPUT PARAMETERS
        if (missing(t)) {
                stop("Enter input parameter t.")
        }

        ## Characteristic function of the (compound) Binomial distribution
        szt <- dim(t)
        t <- c(t)

        if (missing(cfX)) {
                expit <- exp(1i * t)
        } else {
                expit = cfX(t)
        }


        cf <- hypergeom1F1(a,a+b,expit-1)$f

         cf[t==0] <- 1;

         dim(cf) <- szt

        return(cf)
}

