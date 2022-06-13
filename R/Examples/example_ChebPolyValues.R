## EXAMPLE1 (Chebyshev polynomial values specified by coefficients of Sine)
n <- 2^5+1
domain <- c(-pi,pi)
nodes  <- ChebPoints(n,domain)
f <-list(sin(nodes[[1]]))
coeffs <- ChebCoefficients(f)
x <- seq(-pi,pi,length.out=100)
pval <- ChebPolyValues(coeffs,x)
plotReIm(function(x)
           pval,x,
           xlab = "x",
           ylab ="Chebyshev polynomial",
           title ="Chebyshev Polynomial Specified by its Coefficients")

## EXAMPLE2 (Chebyshev polynomial values of the Sine and Cosine functions)
n <- 2^5+1
domain <- c(-pi,pi)
nodes  <- ChebPoints(n,domain)
f <- list(sin(nodes[[1]]), cos(nodes[[1]]), sin(nodes[[1]])*cos(nodes[[1]]))
coeffs <- ChebCoefficients(f);
x <- seq(-pi,pi,length.out=101)
pval_x <- ChebPolyValues(coeffs,c(),domain)
pval_x1<-list()
for (i in 1:length(pval_x[1,])) {
    pval_x1[[i]]   <- pval_x[,i]
  }
plotReIm3(
  pval_x1,x,
  xlab="x",
  ylab="Chebyshev polynomial",
  title="Chebyshev Polynomial Specified by its Coefficients")

## EXAMPLE3 (Chebyshev polynomial values of the Normal PDF and CDF)
n <- 2^7+1
domain <- c(-8, 8)
nodes  <- ChebPoints(n,domain)
f <- list(dnorm(nodes[[1]]), pnorm(nodes[[1]]))
coeffs <- ChebCoefficients(f)
x <- seq(-4,4,length.out=101)
pval_x <- ChebPolyValues(coeffs,x,domain)
pval_x1 <-list()
for (i in 1:length(pval_x[1,])) {
    pval_x1[[i]] <- pval_x[,i]
    }

plotReIm3(
pval_x1,x, title="Chebyshev Polynomial Values Specified by its Coefficients",
xlab = "x",
ylab ="Chebyshev polynomial" )

