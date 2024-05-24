## EXAMPLE 1
a <- 3
b <- 5
c <-
X <- c(1, 2)
f <- Hypergeom1F1MatApprox(a, b, X)

## EXAMPLE 2
# PDF/CDF of minus log Wilks Lambda RV (p=10, n=30, q=5) from its CF
# Here, cf_LogRV_WilksLambdaNC id based on using Hypergeom1F1MatApprox
p <- 10
n <- 30
q <- 5
Delta <- c(1, 2, 3, 10, 50) # nonzero eigenvalues of non-centrality matrix
coef <- -1
cf <- function(t) cf_LogRV_WilksLambdaNC(t, p, n, q, Delta, coef)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, prob = prob, options = options)
