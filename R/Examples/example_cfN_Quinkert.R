## EXAMPLE 1
# CF of the Quinkert distribution with the parameter a=3, b=5
  a <- 3;
  b <- 5;
  t <- seq(from = -15,
           to = 15,
           length.out =501)
  plotReIm(function(t)
          cfN_Quinkert(t, a, b),
          t,
          title = "CF of the Quinkert distribution with a = 3, b = 5")

  ## EXAMPLE 2
  # CF of the compound Quinkert-Exponential distribution
  a <- 3
  b <- 5
  lambda <- 5
  cfX <- function(t)
          cfX_Exponential(t, lambda)
  t <- seq(from = -15,
           to = 15,
           length.out = 501)
  plotReIm(function(t)
          cfN_Quinkert(t, a, b, cfX),
          t,
          title = "CF of the compound Quinkert-Exponential distribution")

 ## EXAMPLE3
 # PDF/CDF of the compound Quinkert-Exponential distribution
    a <- 3;
    b <- 5;
    lambda <- 5;
    cfX <- function(t)
            cfX_Exponential(t, lambda)
    cf <- function(t)
            cfN_Quinkert(t,a,b,cfX)
    x <- seq(0,1.5, length.out = 101)
    prob <- c(0.9, 0.95, 0.99)
    options <- list()
    options$isCompound = TRUE
    result <- cf2DistGP(cf,x,prob,options)




