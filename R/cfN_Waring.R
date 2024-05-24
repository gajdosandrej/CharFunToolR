#' @title Characteristic function of the Waring distribution
#'
#' @description
#' \code{cfN_Waring(t,a,b,r)} evaluates the characteristic function \eqn{cf(t)} of the
#' Waring distribution, with the parameters \eqn{a} (a > 0), \eqn{b} (b > 0), and \eqn{r} (r > 0), i.e.
#'
#'
#' \eqn{cf(t) = cfN_Waring(t,a,b,r)  = ((gamma(a+r)*gamma(a+b)) / (gamma(a)*gamma(a+b+r)))* 2F1(r,b,a+b+r,e^(1i*t))};
#' where 2F1 denotes the Gauss hypergeometric function. The Waring distribution is also known as beta negative binomial distribution. For
#' more details see [4], p. 643
#'
#' @family Discrete Probability Distribution
#'
#' @references
#' [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
#' distribution based on numerical inversion of the compound empirical
#' characteristic function of frequency and severity. Preprint submitted
#' to Insurance: Mathematics and Economics.
#'
#' [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM
#' for computing distributions of collective risk models. Preprint submitted
#' to Journal of Statistical Software.
#'
#' [3] WITKOVSKY V. (2016). Numerical inversion of a characteristic function:
#' An alternative tool to form the probability distribution
#' of output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'
#' [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
#' probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN 3-87773-025-6.
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Beta_negative_binomial_distribution}
#'
#' @param t vector or array of real values, where the CF is evaluated
#' @param a vector of the 'shape' parameters \code{a > 0}. If empty, default value is \code{a = 1}.
#' @param b vector of the 'shape' parameters \code{b > 0}. If empty, default value is \code{b = 1}.
#' @param r number of successes until the experiment is stopped (integer but can be extended to real).
#' @param cfX function.
#'
#' @return Characteristic function \eqn{cf(t)} of the Waring distribution.
#'
#' @note Ver.: 31-Jul-2021 12:47:54 (consistent with Matlab CharFunTool v1.5.1, 15-Nov-2016 13:36:26).
#'
#' @example R/Examples/example_cfN_Waring.R
#'
#' @export
#'
 cfN_Waring <- function(t, a, b, r, cfX) {
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

   cf <- exp(GammaLog(a+r) + GammaLog(a+b) - GammaLog(a) - GammaLog(a+b+r))

   cf <- cf * Hypergeom2F1(r,b,a+b+r,expit)
   cf[t == 0] <- 1

   dim(cf) <- szt

   return(cf)
 }


