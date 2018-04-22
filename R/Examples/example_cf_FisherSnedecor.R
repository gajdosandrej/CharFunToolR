## EXAMPLE 1
# CF of a linear combination of independent F RVs
df1 <- 3:12
df2 <- seq(14, 5, -1)
coef <- 1/10
t <- seq(from = -5, to = 5, length.out = 501)
plotGraf(function(t)
  cf_FisherSnedecor(t, df1, df2, coef), t, title = 'Characteristic function of the linear combination of F RVs')

## EXAMPLE 2
# PDF/CDF  of a linear combination of independent F RVs
df1 <- 3:12
df2 <- seq(14, 5, -1)
coef <- 1/10
cf <- function(t) cf_FisherSnedecor(t, df1, df2, coef)
options <- list()
options$N <- 2^10
options$xMin = 0
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf = cf, prob = prob, options = options)
result
