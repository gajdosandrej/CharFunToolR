#' @title
#'Values of the Chebyshev polynomial
#'
#' @description
#' \code{ChebPolyValues(coeffs,x,domain)} evaluates the values of the Chebyshev polynomial specified
#' by its coefficients at the specified grid of points x, where \eqn{x = xMin +
#' xx*(xMax-xMin)}, with \eqn{domain = [xMin, xMax]} and xx being the grid of
#' points from the interval \eqn{[-1,1]}.
#'
#' @family Utility Function
#'
#' @param coeffs Chebyshev coefficients
#' @param x points, in which we want to compute the value of polynomial.
#' @param domain vector containing lower and upper bound of interval.
#'
#' @return Function returns values of nth order Chebysev polynomial of the first kind evaluated in points \code{x}.
#'
#' @note Ver.: 15-Nov-2021 18:19:46 (consistent with Matlab CharFunTool v1.3.0, 28-May-2021 14:28:24).
#'
#' @example R/Examples/example_ChebPolyValues.R
#'
#' @export

ChebPolyValues <- function(coeffs,x,domain){

        ## CHECK THE INPUT PARAMETERS
        if(missing(domain)) {
                domain <- vector()
        }
        if(missing(x)) {
                x <- vector()
        }
        szx <- dim(x)
        x   <- c(x)

        # Set x and the domain such that x in [-1,1] and domain = [xMin, xMax]
        if(length(domain)==0) {

                if(length(x)==0){
                domain <- c(-1,1)
                x <- seq(-1,1,length.out=101)
                }
                 else{
                domain <- c(min(x),max(x))
                x <- 2*(c(x) - domain[1])/ (domain[2]-domain[1]) - 1
                 }
        }
         else{
                 if (length(x)==0){
                 x <- seq(-1,1,length.out=101)
                 }
                 else{
                         if (min(x) < domain[1] || max(x) > domain[2]){
                                 warning("Some values are outside the stated domain")
                         }
                         x <- 2*(c(x) - domain[1])/ (domain[2]-domain[1]) - 1
                 }

         }
        n<-seq(0,length(coeffs[,1])-1)

        ## ALGORITHM
        s<-matrix(nrow = length(acos(x)),ncol=length(n))

for (i in 1:length(acos(x))) {
        for (j in 1:length(n) ){
           s[i,j]<-   n[j]*acos(x[i])
        }
}

        pval<-cos(s)%*%coeffs

        x<-domain[2]*(x+1)/2+domain[1]*(1-x)/2


      ##  if(all(szx) && dim(pval,2)==1){

        ##        dim(pval) <- szx

        ##}
return(pval)
}


