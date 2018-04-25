## EXAMPLE1
# Calculate CDF/PDF of N(0,1) by inverting its CF
cf <- function(t)
        exp(-t ^ 2 / 2)
result <- cf2DistGP(cf)

## EXAMPLE2
# PDF/CDF of the compound Binomial-Exponential distribution
n <- 25
p <- 0.3
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
x <- seq(0, 5, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE3
# PDF/CDF of the compound Poisson-Exponential distribution
lambda1 <- 10
lambda2 <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda2)
cf <- function(t)
        cfN_Poisson(t, lambda1, cfX)
x <- seq(0, 8, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
