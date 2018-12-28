## EXAMPLE 1
# CF of a linear combination of chi-square random variables
df <- 1
cfX <- function(t) cfX_ChiSquare(t, df)
coef <- 1 / (1:100)
cf <- function(t) cf_Conv(t, cfX, coef)
t <- seq(-10, 10, length.out = 501)
plotReIm(cf, t, title = 'CF of a Linear Combination of iid Chi-Square RVs')
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, options = options)

## EXAMPLE 2
# CF of a linear combination of iid RVs sampled from the empirical
# distribution function (atrtificially generated data)
n <- 30
p <- c(0.2, 0.7, 0.1)
data <- p[1] * rchisq(n, 5) + p[2] * rnorm(n, 10, 1) + p[3] * rt(n, 1)
cfE <- function(t) cfE_Empirical(t, data)
coef <- 1 / c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
cf <- function(t) cf_Conv(t, cfE, coef)
t <- seq(-10, 10, length.out = 501)
plotReIm(cf, t, title = 'CF of a Linear Combination of iid RVs')
result <- cf2DistGP(cf)
