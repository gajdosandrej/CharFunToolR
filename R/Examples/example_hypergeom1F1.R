## EXAMPLE 1
a <- 10
b <- 15
z <- c(0, 1, 10i, 50i, 10 + 50i, 100 + 50i)
n <- 64
result <- hypergeom1F1(a, b, z, n)

## EXAMPLE 2
# CF of Beta(1/2,1/2) distribution
a  = 1 / 2
b  = 1 / 2
t  = seq(-100, 100, length.out = 1001)
plotGraf(
        function(t)
                hypergeom1F1(a, a + b, 1i * t)$f,
        t,
        title = "Characteristic function of Beta(1/2,1/2) distribution",
        xlab = "t",
        ylab = "CF"
)
