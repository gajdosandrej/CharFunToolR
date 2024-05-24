## EXAMPLE 1
# CF of the Gamma distribution with alpha = 2, beta = 2
alpha <- 2
beta <- 2
t <- seq(-20, 20, length.out = 501)
plotReIm(function(t)
        cfX_Gamma(t, alpha, beta),
        t,
        title = "CF of the Gamma distribution with alpha = 2, beta = 2")

## EXAMPLE 2
# PDF/CDF of the Gamma distribution with alpha = 2, beta = 2
alpha <- alpha[18]
beta <- beta[18]
x <- seq(from = 0.1,
         to = 4,
         length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0.1
options$xMax <- 4
options$N <- 2 ^ 14
cf <- function(t)
        cfX_Gamma(t, alpha, beta)
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE 3
# PDF/CDF of the compound Binomial-Gamma distribution
n = 25
p = 0.3
alpha <- 2
beta <- 2
cfX <- function(t)
        cfX_Gamma(t, alpha, beta)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
x <- seq(from = 0,
         to = 25,
         length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
