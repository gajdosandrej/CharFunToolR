## EXAMPLE 1
# CF of the ChiSquared distribution with df = 1
df <- 1
t <- seq(from = -50,to = 50,length.out = 501)
plotGraf(function(t)
  cfX_ChiSquare(t,df), t, title = "CF of the Chi-square distribution with df = 1")

## EXAMPLE 2
# PDF/CDF of the ChiSquare distribution with df = 3
df <- 3
x <- seq(from = 0,to = 15,length.out = 101)
prob <- c(0.9, 0.95, 0.99)
cf <- function(t) cfX_ChiSquare(t,df)
options <- list()
options$xMin <- 0
options$N <- 2^14
result <- cf2DistGP(cf,x,prob,options)

## EXAMPLE 3
# PDF/CDF of the compound Binomial-ChiSquared distribution
n <- 25
p <- 0.3
df <- 3
cfX <- function(t) cfX_ChiSquare(t,df)
cf <- function(t) cfN_Binomial(t,n,p,cfX)
x <- seq(from = 0,to = 80,length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf,x,prob,options)
