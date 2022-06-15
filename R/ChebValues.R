#' @title ChebValues
#'
#' @description \code{ChebValues(coeffs)} converts the n-dimensional Chebyshev coefficients to
# n-dimensional function values, say \eqn{f(x)}, evaluated at the n-dimensional
# Chebyshev points x of the 2nd kind.
#'
#' @family Utility Function
#'
#' @importFrom stats fft
#' @importFrom matlab flipud
#'
#' @param  coeffs Chebyshev coefficients
#'
#' @return  \eqn{V = ChebValues(C)} returns the n-dimensional vector of (polynomial)
#' values evaluated at the Chebyshev points x such that \eqn{V(i) = f(x(i))=
#' C(1)*T_{0}(x(i)) + C(2)*T_{1}(x(i)) + ... + C(N)*T_{N-1}(x(i))} (where \eqn{T_k(x)}
#' denotes the k-th 1st-kind Chebyshev polynomial, and \eqn{x(i)} are the
#' 2nd-kind Chebyshev nodes.
#'
#'If the input C is an (n x m)-matrix then \eqn{V = ChebValues(C)} returns
#'(n x m)-matrix of values V such that \eqn{V(i,j) = P_j(x_i) =
#'C(1,j)*T_{0}(x(i)) + C(2,j)*T_{1}(x(i)) + ... + C(N,j)*T_{N-1}(x(i))}.
#'
#' @note Ver.: 15-Nov-2021 17:55:08 (consistent with Matlab CharFunTool v1.3.0, 28-May-2021 14:28:24).
#'
#' @example R/Examples/example_ChebValues.R
#'
#' @export

ChebValues<-function (coeffs){
    n<-length(coeffs[[1]])
    columns<-ncol(coeffs)
        if ( n <= 1 ){
                values <- coeffs
        }

    isEven<-matrix(data = FALSE,nrow = 1,ncol = columns)
    isOdd  <-matrix(data = FALSE,nrow = 1,ncol = columns)

    abs<-matrix(nrow = n/2,ncol = columns)
    abs1<-matrix(nrow = n/2+1,ncol = columns)
       for (i in 1:columns) {
         abs[,i]<-abs(coeffs[seq(2,n,2),i])
         abs1[,i]<-abs(coeffs[seq(1,n,2),i])
       }

     # Check symmetry
    for (i in 1:columns) {
        if(max(abs[ ,i]) == 0){
                isEven[1,i] <-TRUE
        }
        if(max(abs1[ ,i]) == 0){
                isOdd[1,i]  <-TRUE
        }
    }

        ## Scaling
        coeffs [seq(2,n-1),] <- coeffs[seq(2,n-1),]/2

        ## DCT using an FFT
        tmp<-matrix(nrow = 2*n-2, ncol = columns)

        for (i in 1:columns) {
            tmp[,i]<-c(coeffs[,i],coeffs[seq(n-1,2,-1),i])
        }
        isreal<-TRUE
        isreal1<-TRUE

        for (i in 1: length(coeffs)) {
            if(Im(coeffs[i])!=0){
                isreal<-FALSE
            }
            if(Im(1i*coeffs[i])!=0){
                isreal1<-FALSE
            }
        }
        values<-matrix(nrow = 2*n-2, ncol = columns)


        if (isreal==TRUE){
            for (i in 1: columns) {
                values[,i]<-Re(stats::fft(tmp[,i]))
            }
        }
        else if( isreal==TRUE){              ## Imaginary-valued case
            for (i in 1: columns) {
                values[,i]<-1i*Re(stats::fft(Im(tmp[,i])))
            }
        }

        else{
            for (i in 1: columns) {
                values[,i]<-fft(tmp[,i])
            }
        }


        ## Flip and truncate:
        values1<-matrix(nrow = n,ncol = columns)


        for (i in 1:columns) {
            values1 [ ,i]         <- values[seq(n,1,-1),i]
        }

        ## Symmetry
        values1[,isEven] <- (values1[,isEven]+ matlab::flipud(values1[,isEven]))/2
        values1[,isOdd]  <- (values1[,isOdd] - matlab::flipud(values1[,isOdd]))/2
    return(values1)
        }

