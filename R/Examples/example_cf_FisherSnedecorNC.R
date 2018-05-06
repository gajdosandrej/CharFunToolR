## EXAMPLE 1
# CF of the non-central F RV with delta = 1
df1 <- 3
df2 <- 5
delta <- 1
t <- seq(from = -10,
         to = 10,
         length.out = 201)
plotReIm(function(t)
        cf_FisherSnedecorNC(t, df1, df2, delta), t, title = 'CF of non-central Fisher Snedecor RV')

## EXAMPLE 2
# CDF/PDF of the non-central F RV with delta = 1
df1 <- 3
df2 <- 5
delta <- 1
cf <- function(t)
        cf_FisherSnedecorNC(t, df1, df2, delta)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3
# CDF/PDF of the linear combination of non-central F RVs
df1 <- c(5, 4, 3)
df2 <- c(3, 4, 5)
delta <- c(0, 1, 2)
coef  <- 1 / 3
cf <- function(t)
        cf_FisherSnedecorNC(t, df1, df2, delta, coef)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf = cf, options = options)
