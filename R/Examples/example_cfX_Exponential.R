## EXAMPLE 1
# CF of the Exponential distribution with lambda = 5
lambda <- 5
t <- seq(from = -50,
         to = 50,
         length.out = 501)
plotGraf(function(t)
        cfX_Exponential(t, lambda), t, title = "CF of the Exponential distribution with lambda = 5")

## EXAMPLE 2
# PDF/CDF of the Exponential distribution with lambda = 5
lambda <- 5
cf <- function(t)
        cfX_Exponential(t, lambda)
x  <- seq(from = 0,
          to = 1.5,
          length.out = 101)
options <- list()
options$xMin <- 0
options$SixSigmaRule <- 8
result <- cf2DistGP(cf = cf, x = x, options = options)

## EXAMPLE 3
# PDF/CDF of the compound Binomial-Exponential distribution
n <- 25
p <- 0.3
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
x <- seq(from = 0,
         to = 5,
         length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
