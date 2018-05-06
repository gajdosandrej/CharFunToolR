## EXAMPLE1
# CF of the symmetric Trapezoidal distribution, lambda = 0.5
lambda <- 0.5
t <- seq(-50, 50, length.out = 501)
plotReIm(function(t)
        cfS_Trapezoidal(t, lambda), t,
        title = "CF of the symmetric Trapezoidal distribution with lambda = 0.5")

## EXAMPLE2
# PDF/CDF of the symmetric Trapezoidal distribution, lambda = 0.5
lambda <- 0.5
cf <- function(t)
        cfS_Trapezoidal(t, lambda)
x <- seq(-1, 1, length.out = 100)
xRange <- 2
options <- list()
options$N <- 2 ^ 10
options$dx <- 2 / pi / xRange
result <- cf2DistGP(cf = cf, x = x, options = options)
