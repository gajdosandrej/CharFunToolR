## EXAMPLE1 (CF of the Rectangular distribution on (-2,1))
t <- seq(-50, 50, length.out = 501)
plotReIm(function(t)
        cfX_Rectangular(t, a = -2, b = 1),
        t,
        title = "CF of the Rectangular distribution on (-2,1)")

## EXAMPLE2 (PDF/CDF of the Rectangular distribution on (-2,1))
cf <- function(t)
        cfX_Rectangular(t, a = -2, b = 1)
x <- seq(-2, 1, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
xRange <- 3
options <- list()
options$N <- 2 ^ 10
options$dx <- 2 / pi / xRange
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE3 (PDF/CDF of the Rectangular distribution on (-2,1))
cf <- function(t)
        cfX_Rectangular(t, a = -2, b = 1)
x <- seq(-2, 1, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- -2
options$xMax <- 1
options$N <- 2
result <- cf2DistGP(cf = cf, options = options)
