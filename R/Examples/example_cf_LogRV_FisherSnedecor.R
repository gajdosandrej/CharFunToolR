## EXAMPLE 1
# CF of a weighted linear combination of independent log-F RVs
coef <- c(1, 2, 3, 4, 5)
weight <- coef / sum(coef)
df1 <- 5
df2 <- 3
t <- seq(from = -10,
         to = 10,
         length.out = 201)
plotGraf(function(t)
        cf_LogRV_FisherSnedecor(t, df1, df2, weight),
        t,
        title = "Characteristic function of a linear combination of log-F RVs")

## EXAMPLE 2
# PDF/CDF from the CF by cf2DistGP
coef <- c(1, 2, 3, 4, 5)
weight <- coef / sum(coef)
df1 <- 5
df2 <- 3
cf <- function(t)
        cf_LogRV_FisherSnedecor(t, df1, df2, weight)
options <- list()
options$N <- 2 ^ 12
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf = cf, prob = prob, options = options)
str(result)
