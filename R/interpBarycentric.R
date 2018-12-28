#' @title
#' Barycentric interpolation
#'
#' @description
#' \code{interpBarycentric(x, fun, xNew, options)} evaluates (interpolates) function values
#' \code{funNew} at given points \code{xNew} by barycentric interpolation from function values \code{fun}
#' given at chebpoints \code{x}.
#'
#' @family Utility Function
#'
#' @seealso
#' For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Barycentric_coordinate_system}.
#'
#' @param x points in which is fun given (for more accuracy use chebpoints).
#' @param fun function values of fun given at points \code{x}.
#' @param xNew point in which fun will be estimated.
#' @param options  optional parameter, set the ChebPoints TRUE or FALSE by options$isChebPts <- TRUE.
#'
#' @return This function returns a list consisting of values \code{funNew} of function \code{fun}
#' evaluated at points \code{xNew}.
#'
#' @note Ver.: 16-Sep-2018 21:16:14 (consistent with Matlab CharFunTool v1.3.0, 24-Jul-2017 10:06:48).
#'
#' @example R/Examples/example_interpBarycentric.R
#'
#' @export
interpBarycentric <- function(x, fun, xNew, options) {
  if (missing(xNew)) {
    xNew <- vector()
  }
  if (missing(options)) {
    options <- list()
  }
  if (length(xNew) == 0) {
    xNew <- seq(min(x), max(x))
  }

  if (is.null(options$isChebPts)) {
    options$isChebPts <- TRUE
  }

  x <- c(x)
  fun <- c(fun)

  sztNew <- dim(xNew)
  xNew <- c(xNew)
  nx <- length(x)
  nxNew <- length(xNew)
  funNew <- seq(0, 0, length.out = length(xNew))

  w <- (-1) ^ seq(0, length(x) - 1)
  w[1] <- w[1] / 2
  w[length(x)] <- w[length(x)] / 2

  for (i in seq(1, length(xNew))) {
    A <- 0
    B <- 0
    for (j in seq(1, length(x))) {
      if (xNew[i] == x[j]) {
        exactVal <- TRUE
        funNew[i] <- fun[j]
      } else {
        exactVal <- FALSE
        weight <- w[j] / (xNew[i] - x[j])
        A <- A + weight * fun[j]
        B <- B + weight
      }

      if (exactVal) {
        break
      }
      else {
        funNew[i] <- A / B
      }
    }

  }

  dim(xNew) <- sztNew
  dim(funNew) <- sztNew

  return(list(xNew, funNew))
}
