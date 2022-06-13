## EXAMPLE1 (Chebyshev coefficients of the Sine function on (-pi,pi))
   # n      <- 2^5+1
   # domain <- c(-pi,pi)
   # x      <- ChebPoints(n,domain)
   # f      <-list( sin(x[[1]]))
   # coeffs <- ChebCoefficients(f)
   # print(list(x[[1]], coeffs))
   # x      <- seq(-pi,pi,length.out=100)
   # pval   <- ChebPolyValues(coeffs,x,domain)
   # plotReIm(function(x)
   #         pval,x,
   # xlab="x",
   # ylab="Chebyshev polynomial",
   # title="Chebyshev Polynomial Values Specified by its Coefficients")

## EXAMPLE2 (Chebyshev coefficients of the Sine and the Cosine on (-pi,pi))
   n      <- 2^5+1
   domain <- c(-pi,pi)
   x      <- ChebPoints(n,domain)
   f      <-list(sin(x[[1]]), cos(x[[1]]))
   coeffs <- ChebCoefficients(f)
   print(list(x[[1]], coeffs))
   x      <- seq(-pi,pi,length.out=100)
   pval   <- ChebPolyValues(coeffs,x,domain)
   pval1<-list()
   for (i in 1:length(pval[1,])) {

    pval1[[i]]<-pval[,i]
   }
  # matplot(x, pval,type="o",col=c("red","green"),lty=c(1,1),
   #        xlab="x",
    #       ylab="Chebyshev polynomial",
     #      main="Chebyshev Polynomial Values Specified by its Coefficients")
   plotReIm3( pval1,
              x,
              title="Chebyshev Polynomial Values Specified by its Coefficients",
              xlab = "x",
              ylab ="Chebyshev polynomial" )





