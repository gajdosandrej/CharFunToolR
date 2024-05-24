##EXAMPLE1 (Interpolation by Chebyshev polynomial approximation of Sine)
   n      <- 2^5+1
   domain <- c(-pi,pi)
   nodes  <- ChebPoints(n,domain)
   values <- list (sin(nodes[[1]]))
   x      <- seq(domain[1],domain[2],length.out=100)
   fun    <- InterpChebValues(values,x,domain)
   plotReIm(function(x)
           fun,x,
   xlab="x",
   ylab="Chebyshev Interpolant",
   title="Chebyshev Polynomial Approximation")

## EXAMPLE2 (Interpolation by Chebyshev polynomial approximation)
   n      <- 2^5+1
   domain <- c(-pi,pi)
   nodes  <- ChebPoints(n,domain)
   values <- list(sin(nodes[[1]]), cos(nodes[[1]]), sin(nodes[[1]])*cos(nodes[[1]]))
   x      <- seq(domain[1],domain[2],length.out=100)
   fun<-InterpChebValues(values,x,domain)
   fun1<-list()
   for (i in 1:length(fun[1,])) {
      fun1[[i]]   <- fun[,i]
  }

  plotReIm3(
           fun1,x,
   xlab="x",
   ylab="Chebyshev Interpolant",
   title="Chebyshev Polynomial Approximation")

  ## EXAMPLE3 (CDF/PDF interpolation by Chebyshev polynomial approximation)
     n      <- 2^7+1
     domain <- c(-8,8)
     nodes  <- ChebPoints(n,domain)
     values <- list(dnorm(nodes[[1]]), pnorm(nodes[[1]]))
     x      <- seq(-4,4,length.out=101)
     fun    <- InterpChebValues(values,x,domain)
     fun1<-list()
     for (i in 1:length(fun[1,])) {
        fun1[[i]]   <- fun[,i]
     }
     plotReIm3(fun1,x,
        xlab="x",
        ylab="Chebyshev Interpolant",
        title="Chebyshev Polynomial Approximation")

