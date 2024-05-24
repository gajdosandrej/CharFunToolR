## EXAMPLE 1
# CF of the symmetric Wigner Semicircle distribution on (-1,1)
t <- seq(from = -30,
         to = 30,
         length.out =501)
plotReIm(function(t)
        cf_WignerSemicircle(t),
        t,
        title = "CF of the Wigner Semicircle distribution on (-1,1)")



##EXAMPLE2
# PDF/CDF of Wigner Semicircle  RVs
cf <- function(t)
        cf_WignerSemicircle(t)
x <- seq(-1,1,length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- -1
options$xMax <- 1
result <- cf2DistGP(cf, x, prob, options)

##EXAMPLE3
# CF of a linear combination of independent Wigner Semicircle RVs
t <- seq(-1,1, length.out=  501)
mu <- c(0, 1, 2, 1, 0)
R <- c(1, 1, 2, 2, 3)
coef <- c(1, 2, 3, 4, 5 )
plotReIm(function(t)
        cf_WignerSemicircle(t, mu, R, coef),
        t,
        title = "CF of a linear combination of independent Wigner Semicircle RVs")

## EXAMPLE 4
# PDF/CDF of a linear combination of independent Wigner Semicircle RVs
mu <- c(0, 1, 2, 1, 0)
R <- c(1, 1, 2, 2, 3)
coef <- c(1, 2, 3, 4, 5)
cf <- function(t)
        cf_WignerSemicircle(t, mu, R, coef)
x <- seq(-20, 40, length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N = 2^12
result <- cf2DistGP(cf, x, prob, options)
