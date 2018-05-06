## EXAMPLE 1
# CF of the PearsonV distribution
alpha <- 3 / 2
beta <- 2 / 3
t <- seq(-10, 10, length.out = 1001)
plotReIm(function(t)
        cfX_PearsonV(t, alpha, beta), t,
        title = "CF of the PearsonV distribution with alpha = 3/2, beta = 2/3")

## EXAMPLE 2
# PDF/CDF of the Beta distribution with alpha = 3/2, beta = 2/3
alpha <- 3 / 2
beta <- 2 / 3
prob <- c(0.9, 0.95, 0.99)
x <- seq(0, 40, length.out = 101)
cf <- function(t)
        cfX_PearsonV(t, alpha, beta)
options <- list()
options$xMin <- 0
options$N <- 2 ^ 10
options$SixSigmaRule <- 10
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE 3
# PDF/CDF of the compound Binomial-PearsonV distribution
n <- 25
p <- 0.3
alpha <- 3 / 2
beta <- 2 / 3
prob <- c(0.9, 0.95, 0.99)
x <- seq(0, 200, length.out = 101)
cfX <- function(t)
        cfX_PearsonV(t, alpha, beta)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
options <- list()
options$isCompound <- TRUE
options$N <- 2 ^ 10
result <- cf2DistGP(cf, x, prob, options)
