## EXAMPLE 1
# CF of the symmetric Triangular distribution on (-1,1)
t <- seq(from = -50,
         to = 50,
         length.out = 501)
plotReIm(function(t)
        cfS_Triangular(t), t, title = "CF of the the symmetric Triangular distribution on (-1,1)")

## EXAMPLE 2
# PDF/CDF of the the symmetric Triangular distribution on (-1,1)
cf <- function(t)
        cfS_Triangular(t)
x <- seq(from = -1,
         to = 1,
         length.out = 101)
xRange <- 2
options <- list()
options$dt <- 2 * pi / xRange
result <- cf2DistGP(cf = cf, x = x, options = options)
