## EXAMPLE 1
# CF of the log-transformed non-central F RV with delta = 1 and coef = -1
df1 <- 3
df2 <- 5
delta <- 1
coef  <- -1
t <- seq(from = -10, to = 10, length.out = 201)
plotGraf(function(t)
  cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef), t, title = 'CF of minus log-transformed F RV')

## EXAMPLE 2
# CDF/PDF of the minus log-transformed non-central F RV with delta = 1
df1 <- 3
df2 <- 5
delta <- 1
coef <- -1
cf <- function(t) cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef)
options <- list()
options$N  <- 2^12
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3
# CDF/PDF of the linear combination of log-transformed non-central F RVs
df1 <- c(5, 4, 3)
df2 <- c(3, 4, 5)
delta <- c(0, 1, 2)
coef  <- -1/3
cf <- function(t) cf_LogRV_FisherSnedecorNC(t, df1, df2, delta, coef)
options <- list()
options$N  <- 2^12
result <- cf2DistGP(cf = cf, options = options)
