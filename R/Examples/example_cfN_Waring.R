## EXAMPLE1
# CF of the Waring distribution with a = 2.2, b = 3.3, r = 4
# The CF is not computed correctly!! Because the hypergeomtric function
# Hypergeom2F1(r,b,a+b+r,z) does not converege abs(z)>=1. Here z = exp(1i*t),
# ans abs(exp(1i*t)) = 1.
 a <- 2.2
 b <- 3.3
 r <- 4
 t <- seq(from = -5,
               to = 5,
               length.out =1001)
 plotReIm(function(t)
         cfN_Waring(t, a, b, r),
         t,
         title = "CF of the Waring distribution with a = 2.2, b = 3.3, r = 4")



##EXAMPLE2
 #CF of the compound Waring-Exponential distribution
 a <- 2.2
 b <- 3.3
 r <- 4
 lambda <- 5
 cfX <- function(t)
         cfX_Exponential(t,lambda)
 t <- seq(-10,10,length.out = 501)
 plotReIm(function(t)
          cfN_Waring(t, a, b, r, cfX),
          t,
           title = "CF of the compound Waring-Exponential distribution")

##EXAMPLE3
# PDF/CDF of the compound Waring-Exponential distribution
  a <- 2.2
  b <- 3.3
  r <- 4
  lambda <- 5;
  cfX <- function(t)
          cfX_Exponential(t,lambda)
  cf <- function(t)
          cfN_Waring(t, a, b, r, cfX)
  x <- seq(0,35, length.out=  101)
  prob <- c(0.9, 0.95, 0.99)
  options <- list()
  options$isCompound = TRUE
  result <- cf2DistGP(cf,x,prob,options)
