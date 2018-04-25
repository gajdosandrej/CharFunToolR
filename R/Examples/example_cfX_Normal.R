## EXAMPLE1 (CF of the Normal distribution N(1,1))
t <- seq(-5, 5, length.out = 501)
plotGraf(function(t)
        cfX_Normal(t, mean = 1, variance = 1), t, title = "CF of the Normal distribution N(1,1)")

## EXAMPLE2 (PDF/CDF of the Normal distribution N(1,1))
cf <- function(t)
        cfX_Normal(t, mean = 1, variance = 1)
x <- seq(-4, 4, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2 ^ 5
options$SixSigmaRule <- 8
result <- cf2DistGP(cf, x, prob, options)
