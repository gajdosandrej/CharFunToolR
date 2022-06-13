## EXAMPLE1
# CF of the Laplace RV
mu <- 0
beta <- 1
t  <- seq(from = -10,
              to = 10,
              length.out =201)
plotReIm(function(t)
        cf_Laplace(t, mu, beta),
        t,
        title = "Characteristic function of the Laplace RVs")



##EXAMPLE2
# PDF/CDF of the Laplace RV
mu <- 0
beta <- 1
x <- seq(-5, 5, length.out = 101)

prop <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf <- function(t)
        cf_Laplace(t, mu, beta)
result <- cf2DistGP(cf, x, prob)

##EXAMPLE3
# PDF/CDF of the linear combination of Laplace RVs
mu <- c(-4, -1, 2, 3)
beta <- c(0.1, 0.2, 0.3, 0.4)
coef <- c(1, 2, 3, 4)
prob <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf <- function(t)
        cf_Laplace(t, mu, beta, coef)
options <- list()
options$N <- 2^12
result <- cf2DistGP(cf,prob=prob,options=options)

## EXAMPLE 4
# PDF/CDF of the linear combination of the Laplace RVs
mu <- c(-10, 10, 20, 30, 40)
beta <- c(1, 2, 3, 4, 5)
coef <- c(1/2, 1, 3/4, 5, 1)
prop <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf <- function(t)
        cf_Laplace(t, mu, beta, coef)
options <- list()
options$N <- 2^12
result <- cf2DistGP(cf,c(),prob,options)
