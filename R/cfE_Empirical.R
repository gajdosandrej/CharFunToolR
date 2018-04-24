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

  ##
  if (length(data) == 0) {
    data <- 1
  }

  weights <- 1 / length(data)
  data <- c(data)

  # Special treatment for mixtures with large number of variables

  szcoefs <- dim(t(data))
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

    if (missing(cfX)) {
      aux <- exp(1i * t %*% t(data[idx]))
    } else {
      aux <- t(apply(as.matrix(cfX(t)), 1, FUN = '^', data[idx]))
    }

    if (length(weights) == 1) {
      cf <- cf + apply(weights * aux, 1, sum)
    } else {
      cf <- cf + apply(aux * weights[idx], 1, sum)
      cf <- cf + apply(t(apply(aux * weights[idx], 1, FUN = '*', c(1,2,3))), 1, sum)
    }

  }

  dim(cf) <- szt

  return(cf)
}
