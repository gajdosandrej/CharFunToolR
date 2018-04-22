## EXAMPLE 1
# CF of a weighted linear combination of independent log-ChiSquare RVs
coef <- c(1, 2, 3, 4, 5)
weight <- coef / sum(coef)
df <- c(1, 2, 3, 4, 5)
t <- seq(from = -20, to = 20, length.out = 1001)
plotGraf(function(t)
  cf_LogRV_ChiSquare(t, df, weight), t, title = "CF of a linear combination of minus log-ChiSquare RVs")

## EXAMPLE 2
# PDF/CDF of a linear combination of independent log-ChiSquare RVs
coef <- c(1, 2, 3, 4, 5)
weight <- coef / sum(coef)
df <- c(1, 2, 3, 4, 5)
cf <- function(t) cf_LogRV_ChiSquare(t, df, weight)
options <- list()
options$N <- 2^12
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf = cf, prob = prob, options = options)
result
