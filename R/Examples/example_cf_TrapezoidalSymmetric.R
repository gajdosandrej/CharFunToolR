## EXAMPLE 1
#CF of the symmetric Trapezoidal distribution with lambda = 1/2
lambda <- 1 / 2
t <- seq(from = -50,
         to = 50,
         length.out = 201)
plotReIm(function(t)
        cf_TrapezoidalSymmetric(t, lambda), t,
        title = "CF of the symmetric Trapezoidal distribution on (-1,1)")

## EXAMPLE 2
# CF of a linear combination of independent Trapezoidal RVs
t <- seq(from = -20,
         to = 20,
         length.out = 201)
lambda <- c(3, 3, 4, 4, 5) / 7
coef <- c(1, 2, 3, 4, 5) / 15
plotReIm(function(t)
        cf_TrapezoidalSymmetric(t, lambda, coef), t,
        title = "CF of a linear combination of independent Trapezoidal RVs")

## EXAMPLE 3
# PDF/CDF of a weighted linear combination of independent Trapezoidal RVs
lambda <- c(3, 3, 4, 4, 5) / 7
coef <- c(1, 2, 3, 4, 5) / 15
cf <- function(t)
        cf_TrapezoidalSymmetric(t, lambda, coef)
x <- seq(from = -1,
         to = 1,
         length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2 ^ 12
options$xMin <- -1
options$xMax <- 1
result <- cf2DistGP(cf, x, prob, options)
