#' @title ChebCoefficients
#'
#' @description
#' \code{ChebCoefficients(values)} converts function values, say values = \eqn{f(x)}, evaluated
#' at Chebyshev points of the 2nd kind, say n-dimensional x, to
#' n-dimensional Chebyshev coefficients.
#'
#' @family Utility Function
#'
#' @importFrom pracma ifft
#' @importFrom matlab flipud
#'
#' @param values function values evaluated at Chebyshev points of the 2nd kind
#'
#' @return \eqn{c = ChebCoefficients(values)} returns the n-dimensional vector of coefficients
#' such that \eqn{f(x) = C(1)*T_{0}(x) + C(2)*T_{1}(x)+ ... +
#' C(N)*T_{n-1}(x)} (where \eqn{T_{k}(x)} denotes the k-th 1st-kind Chebyshev
#' polynomial) interpolates the data \eqn{[values(1); ... ; values(n)]} at Chebyshev
#' points \eqn{x = [x(1); ... ; x(n)]} of the 2nd kind.
#'
#' If values is an (n x m)-matrix, then \eqn{C = ChebCoefficients(values)} returns the
#' (n x m)-matrix of coefficients C.
#'
#' @note Ver.: 15-Nov-2021 17:33:19 (consistent with Matlab CharFunTool v1.3.0, 28-May-2021 14:28:24).
#'
#' @example R/Examples/example_ChebCoefficients.R
#'
#' @export


ChebCoefficients <- function(values){

        n<-length(values[[1]])
        columns<-length(values)
        values1<-matrix(nrow =n,ncol = columns )
        for (i in 1:columns) {
          values1[ ,i]<-values[[i]]
        }
        if ( n <= 1 ){

                coeffs <- values
        }
isEven<-matrix(data = FALSE,nrow = 1,ncol = columns)
isOdd  <-matrix(data = FALSE,nrow = 1,ncol = columns)

# Check symmetry
for (i in 1:columns) {


         if(max(abs(values1-matlab::flipud(values1))[ ,i]) ==0){

                isEven[1,i] <-TRUE
        }
        if(max(abs(values1+matlab::flipud(values1))[ ,i]) ==0){
        isOdd[1,i]  <-TRUE
       }
}

tmp<-matrix(nrow = 2*n-2, ncol = columns)
        for (i in 1:columns) {
                    tmp[,i]<-c(values1[seq(n,2,-1),i],values1[seq(1,n-1),i])
        }

        isreal<-TRUE
        isreal1<-TRUE

        for (i in 1: length(values1)) {
        if(Im(values1[i])!=0){
        isreal<-FALSE
        }
                if(Im(1i*values1[i])!=0){
                        isreal1<-FALSE
                }
        }

        coeffs<-matrix(nrow = 2*n-2, ncol = columns)


        if (isreal==TRUE) {
                for (i in 1: columns) {
                        coeffs[,i]<-pracma::ifft(tmp[,i])
                }

                coeffs <- Re(coeffs)
        }
       else if(isreal1==TRUE ){
               for (i in 1: columns) {
                       coeffs[,i]<-pracma::ifft(Im(tmp)[,i])
               }

                coeffs <- 1i*Re(coeffs)
        }
        else{
                for (i in 1: columns) {
                        coeffs[,i]<-pracma::ifft(tmp[,i])
                }
        }

        # Truncate, scale the interior coefficients, and adjust for symmetry
coeffs1<-matrix(nrow = n,ncol = columns)
for (i in 1:columns) {
     coeffs1  [ ,i]         <- coeffs[seq(1,n),i]
}


        coeffs1[seq(2,n-1),]    <- 2*coeffs1[seq(2,n-1),]
        coeffs1[seq(2,n,2),isEven] <- 0
        coeffs1[seq(1,n,2),isOdd]    <- 0
        return(coeffs1)
        }

