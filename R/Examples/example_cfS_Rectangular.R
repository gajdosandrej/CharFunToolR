## EXAMPLE1
# CF of the Rectangular distribution on (-1,1)
t <- seq(-50, 50, length.out = 501)
plotReIm(function(t)
        cfS_Rectangular(t), t, title = "CF of the Rectangular distribution on (-1,1)")

## EXAMPLE2
# PDF/CDF of the Rectangular distribution on (-1,1)
cf <- function(t)
        cfS_Rectangular(t)
x <- seq(-2, 2, length.out = 101)
xRange <- 2
options <- list()
options$N <- 2 ^ 5
options$dx <- 2 / pi / xRange
result <- cf2DistGP(cf = cf, x = x, options = options)
