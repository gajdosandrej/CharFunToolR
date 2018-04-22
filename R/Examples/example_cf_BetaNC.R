## EXAMPLE 1
# CF of non-central Beta RV with delta = 1
alpha <- 1
beta  <- 3
delta <- 1
t <- seq(from = -50, to = 50, length.out = 201)
plotGraf(function(t)
  cf_BetaNC(t, alpha, beta, delta), t, title = 'CF of a Beta RVs')

## EXAMPLE 2
# CDF/PDF of non-central Beta RV with delta = 1
alpha <- 1
beta  <- 3
delta <- 1
cf <- function(t) cf_BetaNC(t, alpha, beta, delta)
options <- list()
options$xMin <- 0
options$xMax <- 1
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3
# CDF/PDF of the linear combination of non-central Beta RVs
alpha <- c(5, 4, 3)
beta  <- c(3, 4, 5)
delta <- c(0, 1, 2)
coef  <- 1/3
cf <- function(t) cf_BetaNC(t, alpha, beta, delta, coef)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf = cf, options = options)
