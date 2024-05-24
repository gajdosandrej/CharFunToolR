## EXAMPLE 1:
# CF of the Logistic RV
mu   <- 0
beta <- 1
t    <- seq(-10, 10, length.out=201)
cf   <- cf_Logistic(t, mu, beta)
plotReIm(function(t) cf_Logistic(t, mu, beta), t,
title = 'Characteristic function of the Logistic RVs')

## EXAMPLE 2:
# CDF of the Logistic RV
mu  <- 0
beta <- 1
x   <- seq(-10,10,length.out=101)
prob <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf  <- function(t) cf_Logistic(t,mu,beta)
result <- cf2DistGP(cf,x,prob)

## EXAMPLE 3:
# PDF/CDF of the linear combination of Logistic RVs
mu   <- c(-4, -1, 2, 3)
beta <- c(0.1, 0.2, 0.3, 0.4)
coef <- c(1, 2, 3, 4)
prob <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf   <- function(t) cf_Logistic(t,mu,beta,coef)
options <- list()
options$N <- 2^12
result <- cf2DistGP(cf,c(),prob,options)

## EXAMPLE 4:
# PDF/CDF of the linear combination of Logistic RVs
mu   <- c(-10, 10, 20, 30, 40)
beta <- c(1, 2, 3, 4, 5)
coef <- c(1/2, 1, 3/4, 5, 1)
prob <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf   <- function(t) cf_Logistic(t,mu,beta,coef)
options <- list()
options$N <- 2^12
result <- cf2DistGP(cf,c(),prob,options)


