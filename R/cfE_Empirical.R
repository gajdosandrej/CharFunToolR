#' @title Characteristic function of the EMPIRICAL distribution, based on the observed data
#'
#' @description
#' That is, \eqn{cf(t) = (1/N)*(cfD(data_1*t) +...+ cfD(data_N*t))}, where
#' \eqn{cfD(t)} represents the characteristic function of the DIRAC RV
#' concentrated at the constant \eqn{d=1}, i.e. \eqn{cfD(t) = exp(1i*t)}.
#'
#' \code{cfE_Empirical(t, data, cfX)} evaluates the compound characteristic function
#' \eqn{cf(t) = cfE_Empirical(-1i*log(cfX(t)),data) = (1/N) * sum_{j=1}^N cfX(t)^data_j};
#' where \code{cfX} is function handle of the characteristic function \eqn{cfX(t)} of the random variable \eqn{X}.
#'
#' @family Empirical Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Empirical_characteristic_function}.
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param data vector of data, i.e. constants where the DIRAC RVs
#' are concentrated. If empty, default value is \code{data = 1}.
#' @param cfX function handle of the characteristic function of a random
#' variable \eqn{X}. If \code{cfX} is non-empty, a compound CF is evaluated as \eqn{cf(t) = cf(-1i*log(cfX(t)),data)}.
#'
#' @return Characteristic function \eqn{cf(t)} of the EMPIRICAL distribution, based on the observed data.
#'
#' @references
#' WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
#' loss distribution based on numerical inversion of the compound empirical
#' characteristic function of frequency and severity. arXiv preprint arXiv:1701.08299.
#'
#' @note Ver.: 23-Sep-2018 14:40:43 (consistent with Matlab CharFunTool v1.3.0, 15-Sep-2018 12:56:00).
#'
#' @example R/Examples/example_cfE_Empirical.R
#'
#' @export
#'
cfE_Empirical <- function(t, data, cfX) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(t)) {
    stop("Enter input parameter t.")
  }
  if (missing(data)) {
    data <- vector()
  }
  if(missing(cfX)) {
          cfX <- NULL
  }

  weights <- vector()
  cf <- cfE_DiracMixture(t, data, weights, cfX)

  return(cf)
}
