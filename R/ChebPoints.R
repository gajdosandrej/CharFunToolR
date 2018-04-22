#' @title
#' Chebpoints on the given interval
#'
#' @description
#' \code{ChebPoints(N, interval)} evaluates Evaluate \eqn{N} chebpoints on the given interval \eqn{[a,b]}.
#'
#' @family
#'
#' @param N number of chebpoints required.
#' @param interval vector containing lower and upper bound of interval.
#'
#' @return  Function returns \eqn{N} chebpoints on the given interval.
#'
#' @example R/Examples/example_ChebPoints.R
#'
#' @export
#'
ChebPoints <- function(N, interval) {
  ## CHECK THE INPUT PARAMETERS
  if (missing(interval)) {
    interval <- vector()
  }
  if (length(interval) == 0) {
    interval <- c(-1, 1)
  } else if (length(interval) != 2) {
    stop(
      "Input parameter interval should be a vector of length 2, containing lower and upper bound."
    )
  }

  a <- interval[1]
  b <- interval[2]
  pts <- (a + b) / 2 - cos(pi * (0:N) / N) * (b - a) / 2

  return(pts)
}
