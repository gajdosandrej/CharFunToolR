## EXAMPLE 1
t <- seq(-10, 10, length.out = 101)
z <- 1i * t
p <- 10
f <- GammaMultiLog(z, p)

## EXAMPLE 2
t <- seq(-10, 10, length.out = 101)
z <- 1i * t
p <- 10
funmode <- 1
f <- GammaMultiLog(z, p, funmode)
