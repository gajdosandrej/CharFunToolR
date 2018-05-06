## EXAMPLE 1
# CF of the minus log-transformed non-central ChiSquare RV with delta = 1
df <- 3
delta <- 1
coef  <- -1
t <- seq(from = -10,
         to = 10,
         length.out = 201)
plotReIm(function(t)
        cf_LogRV_ChiSquareNC(t, df, delta, coef),
        t,
        title = 'CF of minus log-transformed ChiSquare RV')

## EXAMPLE 2
# CDF/PDF of the minus log-transformed non-central ChiSquare RV
df <- 3
delta <- 1
coef <- -1
cf <- function(t)
        cf_LogRV_ChiSquareNC(t, df, delta, coef)
options <- list()
options$N <- 2 ^ 10
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3
# CDF/PDF of the linear combination of the minus log-transformed non-central ChiSquare RVs
df <- c(3, 4, 5)
delta <- c(0, 1, 2)
coef  <- c(-1,-1,-1) / 3
cf <- function(t)
        cf_LogRV_ChiSquareNC(t, df, delta, coef)
options <- list()
options$N <- 2 ^ 10
result <- cf2DistGP(cf = cf, options = options)
