## EXAMPLE1 (CF of the Triangular distribution on (-1 , 2) mode = 0)
t <- seq(-10, 10, length.out = 501)
plotGraf(function(t)
  cfX_Triangular(t, a = -1, b = 2, c = 0), t,
  title = "CF of the Triangular distribution on (-1 , 2)")

## EXAMPLE2 (PDF/CDF of the Triangular distribution on (-1 , 2) mode = 0)
cf <- function(t)
  cfX_Triangular(t, a = -1, b = 2, c = 0)
x <- seq(-1, 2, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
xRange <- 3
options <- list()
options$N <- 2 ^ 10
options$dx <- 2 / pi / xRange
result <- cf2DistGP(cf = cf, x = x, options = options)
