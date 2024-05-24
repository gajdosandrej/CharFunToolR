#' @title Chebpoints on the given interval
#'
#' @description
#' \code{ChebPoints(N, interval)} evaluates \eqn{N} chebpoints on the given interval \eqn{[a,b]}.
#'
#' @family Utility Function
#'
#' @importFrom pracma ifft
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
  }

  if (N == 0){
    x<-vector()
    w <- vector()
  }
  else if( N == 1){
    x <- 0
    w <- 2
  }
  else{
    m<- N - 1
    x <- sin(pi*(seq(-m, m, 2))/(2*m))
    v <- 2/t(c(1,1-(seq(2,N-1,2))^2))
    v <- c(v,v[seq(floor(N/2),2,-1)])
    w <- pracma::ifft(v)
    w[1] <- w[1]/2
    w[N]<-w[1]
  }

  if (interval[1] != -1 || interval[2] != 1){
    x <- interval[2]*(x + 1)/2 + interval[1]*(1 - x)/2
    w <- (diff(interval)/2)*w
  }



  pts <-list(x,w)

  return(pts)
}
