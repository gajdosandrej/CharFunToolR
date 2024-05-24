## EXAMPLE1
# CF of the Polya-Eggenberger distribution with a = 2.2, b = 3.3, m = 4
a <- 2.2
b <- 3.3
m <- 4
t <- seq(-15, 15, length.out = 1001)
plotReIm(function(t)
        cfN_PolyaEggenberger(t, a, b, m), t,
        title = "CF of the Polya-Eggenberger distribution with a = 2.2, b = 3.3, m = 4")

## EXAMPLE2
# CF of the compound Polya-Eggenberger-Exponential distribution
a <- 2.2
b <- 3.3
m <- 4
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
t <- seq(-50, 50, length.out = 501)
plotReIm(function(t)
        cfN_PolyaEggenberger(t, a, b, m, cfX), t,
        title = "CF of the compound Polya-Eggenberger-Exponential distribution")

## EXAMPLE3
# PDF/CDF of the compound Polya-Eggenberger-Exponential distribution
a <- 2.2
b <- 3.3
m <- 4
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
cf <- function(t)
        cfN_PolyaEggenberger(t, a, b, m, cfX)
x <- seq(0, 2.5, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
