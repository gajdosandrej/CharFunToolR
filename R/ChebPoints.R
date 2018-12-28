#' @title
#' Chebpoints on the given interval
#'
#' @description
#' \code{ChebPoints(N, interval)} evaluates \eqn{N} chebpoints on the given interval \eqn{[a,b]}.
#'
#' @family Utility Function
#'
#' @param N number of chebpoints required.
#' @param interval vector containing lower and upper bound of interval.
#'
#' @return  Function returns \eqn{N} chebpoints on the given interval.
#'
#' @note Ver.: 16-Sep-2018 21:16:55 (consistent with Matlab CharFunTool v1.3.0, 24-Jul-2017 10:06:48).
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
