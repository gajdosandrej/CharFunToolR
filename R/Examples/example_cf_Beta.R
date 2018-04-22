## EXAMPLE 1
# CF of a Beta RV
alpha <- 1
beta  <- 3
t <- seq(from = -50, to = 50, length.out = 501)
plotGraf(function(t)
  cf_Beta(t, alpha, beta), t, title = "CF of a Beta RVs")

## EXAMPLE 2
# PDF/CDF of a Beta RV
alpha <- 1
beta  <- 3
cf <- function(t) {cf_Beta(t, alpha, beta)}
options <- list()
options$N <- 2^12
options$xMin <- 0
options$xMax <- 1
x <- seq(from = 0, to = 1, length.out = 201)
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE 3
# CF of a linear combination of independent Beta RVs
alpha <- 1
beta  <- 3
coef  <- 1/(((1:50) - 0.5)*pi)^2
weights <- coef/sum(coef)
t <- seq(from = -100, to = 100, length.out = 501)
plotGraf(function(t)
  cf_Beta(t,alpha,beta,weights), t, title = "CF of a weighted linear combination of independent Beta RVs")

## EXAMPLE 4
# PDF/CDF of a weighted linear combination of independent Beta RVs
alpha <- 1
beta  <- 3
coef  <- 1/(((1:50) - 0.5)*pi)^2
weights <- coef/sum(coef)
cf <- function(t) {cf_Beta(t,alpha,beta,weights)}
options <- list()
options$xMin <- 0
options$xMax <- 1
x <- seq(from = 0, to = 1, length.out = 201)
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf,x,prob,options)
