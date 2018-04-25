## EXAMPLE 1
# CF of the symmetric Beta distribution with theta = 3/2 on (-1,1)
theta <- 3 / 2
t <- seq(from = -50,
         to = 50,
         length.out = 201)
plotGraf(function(t)
        cf_BetaSymmetric(t, theta),
        t,
        title = "CF of the symmetric Beta distribution on (-1,1)")

## EXAMPLE 2
# CF of a linear combination of independent Beta RVs
t <- seq(from = -20,
         to = 20,
         length.out = 201)
theta <- c(3, 3, 4, 4, 5) / 2
coef <- c(1, 2, 3, 4, 5) / 15
plotGraf(function(t)
        cf_BetaSymmetric(t, theta, coef),
        t,
        title = "CF of a linear combination of independent Beta RVs")

## EXAMPLE 3
# PDF/CDF of a weighted linear combination of independent Beta RVs
theta <- c(3, 3, 4, 4, 5) / 2
coef <- c(1, 2, 3, 4, 5) / 15
cf   <- function(t)
        cf_BetaSymmetric(t, theta, coef)
x <- seq(from = -1,
         to = 1,
         length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2 ^ 12
options$xMin <- -1
options$xMax <- 1
result <- cf2DistGP(cf, x, prob, options)
