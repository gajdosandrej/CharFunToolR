## EXAMPLE1 (
# CF of the symmetric Arcsine distribution on (-1,1)
t <- seq(-50, 50, length.out = 501)
plotGraf(function(t)
  cfS_Arcsine(t), t, title = "CF of the the Arcsine distribution on (-1,1)")

## EXAMPLE2
# PDF/CDF of the symmetric Arcsine distribution on (-1,1)
cf <- function(t)
  cfS_Arcsine(t)
x <- seq(-1, 1, length.out = 501)
xRange <- 2
options <- list()
options$N <- 2 ^ 12
options$dt <- 2 * pi / xRange
result <- cf2DistGP(cf = cf, x = x, options = options)
