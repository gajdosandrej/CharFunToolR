## EXAMPLE1
# CF of the symmetric Beta distribution with theta = 3/2 on (-1,1)
theta <- 3 / 2
t <- seq(-50, 50, length.out = 501)
plotGraf(function(t)
  cfS_Beta(t, theta), t, title = "CF of the symmetric Beta distribution on (-1,1)")

## EXAMPLE2
# PDF/CDF of the the symmetric Beta distribution on (-1,1)
theta <- 3 / 2
cf <- function(t)
  cfS_Beta(t, theta)
x <- seq(-1, 1, length.out = 101)
xRange <- 2
options <- list()
options$dx <- 2 * pi / xRange
options$N <- 2 ^ 8
result <- cf2DistGP(cf = cf, x = x, options = options)
