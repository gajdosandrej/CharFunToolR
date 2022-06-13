## EXAMPLE 1
# CF of the linear combination of symmetric Laplace RV
beta <- 1
t <- seq(from = -10,
         to = 10,
         length.out =201)
plotReIm(function(t)
        cfS_Laplace(t, beta),
        t,
        title = "Characteristic function of the symmetric Laplace RV")



##EXAMPLE2
# PDF/CDF of the symmetric Laplace RVs
beta <- 1
x <- seq(-5, 5, length.out = 101)
prop <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf <- function(t)
        cfS_Laplace(t, beta)
result <- cf2DistGP(cf, x, prob)

##EXAMPLE3
# PDF/CDF of the linear combination of symmetric Laplace RVs
beta <- c(0.1, 0.2, 0.3, 0.4)
coef <- c(1, 2, 3, 4)
prop <- c(0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999)
cf <- function(t)
        cfS_Laplace(t, beta, coef)
options <- list()
options$N <- 2^12
result <- cf2DistGP(cf, c(), prob, options)
