## EXAMPLE1
# CF of the Lognormal distribution with mu = 0,sigma = 1
mu <- 0
sigma <- 1
t <- seq(-20, 20, length.out = 2 ^ 10 + 1)
plotGraf(function(t)
  cfX_LogNormal(t, mu, sigma), t, title = "Characteristic function of the Lognormal distribution")

## EXAMPLE2
# CDF/PDF of the Lognormal distribution with mu = 0,sigma = 1
mu <- 0
sigma <- 1
x <- seq(0, 15, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
cf <- function(t) cfX_LogNormal(t, mu, sigma)
options <- list()
options$xMin <- 0
options$N <- 2^13
options$SixSigmaRule <- 8
result <- cf2DistGP(cf,x,prob,options)

## EXAMPLE3
# PDF/CDF of the compound Poisson-Lognormal distribution
mu <- 0
sigma <- 1
lambda <- 10
x <- seq(0, 70, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
cfX <- function(t) cfX_LogNormal(t, mu, sigma)
cf <- function(t) cfN_Poisson(t, lambda, cfX)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
