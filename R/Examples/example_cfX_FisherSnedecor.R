## EXAMPLE 1
# CF of the F-distribution with df1 = 3, df2 = 5
df1 <- 3
df2 <- 5
t <- seq(-30, 30, length.out = 2^10+1)
plotReIm(function(t) cfX_FisherSnedecor(t, df1, df2), t,
         title = 'Characteristic function of the F-distribution')

## EXAMPLE 2
# PDF/CDF of the F-distribution with df1 = 3, df2 = 5
df1 <- 3
df2 <- 5
x <- seq(0, 25, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
cf <- function(t) cfX_FisherSnedecor(t, df1, df2)
options <- list()
options$xMin <- 0
options$xMax <- 500
options$N <- 2^15
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE 3
# PDF/CDF of the compound Binomial-Fisher-Snedecor distribution
n <- 25
p <- 0.3
df1 <- 3
df2 <- 5
cfX <- function(t) cfX_FisherSnedecor(t, df1, df2)
cf <- function(t) cfN_Binomial(t, n, p, cfX)
x <- seq(0, 80, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
str(result)
