## EXAMPLE 1
# CF of the symmetric Arcsine distribution on (-1,1)
t <- seq(from = -50,
         to = 50,
         length.out = 501)
plotReIm(function(t)
        cf_ArcsineSymmetric(t),
        t,
        title = "CF of the the Arcsine distribution on (-1,1)")

## EXAMPLE 2
# CF of a linear combination of independent Arcsine RVs
t <- seq(from = -1,
         to = 1,
         length.out = 501)
coef <- c(1, 2, 3, 4, 5)
plotReIm(function(t)
        cf_ArcsineSymmetric(t, coef),
        t,
        title = "CF of a linear combination of independent Arcsine RVs")

## EXAMPLE 3
## PDF/CDF of a linear combination of independent Arcsine RVs
coef <- c(1, 2, 3, 4, 5)
cf   <- function(t)
        cf_ArcsineSymmetric(t, coef)
x    <- seq(from = -20,
            to = 20,
            length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2 ^ 12
result <- cf2DistGP(cf, x, prob, options)
