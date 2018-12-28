## EXAMPLE 1
# CF of log of noncentral Wilks Lambda RV distribution Lambda(p,n,q,delta)
p <- 5
n <- 10 # d.f. of within SS&P
q <- 3  # d.f. of between SS&P
delta <- sort(runif(p))
coef <- 1
t <- seq(-10, 10, length.out = 201)
# MAX <- 0
cf1 <- function(t) cf_LogRV_WilksLambdaNC(t, p, n, q, delta, coef, 0)
# MAX <- 20
cf2 <- function(t) cf_LogRV_WilksLambdaNC(t, p, n, q, delta, coef, 20)
plotReIm2(list(cf1, cf2), list(t, t), title = 'CF of log of non-central Wilks Lambda RV')

## EXAMPLE 2
# CF of a weighted linear combination of minus log Wilks Lambda RVs
p <- c(5, 5, 5)
n <- c(10, 15, 20)
q <- c(3, 2, 1)
# delta <- c(sort(runif(p[1])), sort(runif(p[2])), sort(runif(p[3])))
delta <- matrix(list(runif(p[1]), runif(p[2]), runif(p[3])), 1, 3)
coef <- -c(10, 15, 20) / 45
t <- seq(-20, 20, length.out = 201)
plotReIm(function(t) {cf_LogRV_WilksLambdaNC(t, p, n, q, delta, coef)}, t,
         title = 'CF of a weighted linear combination of -log Wilks Lambda RVs')

## EXAMPLE 3
# PDF/CDF of minus log Wilks Lambda RV (p=10, n=30, q=5) from its CF
p <- 10
n <- 30
q <- 5
delta <- c(1, 2, 3, 10, 30)
coef <- -1
cf0 <- function(t) cf_LogRV_WilksLambdaNC(t, p, n, q, coef = coef)
cf <- function(t) cf_LogRV_WilksLambdaNC(t, p, n, q, delta, coef)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result0 <- cf2DistGP(cf0, prob = prob, options = options)
result <- cf2DistGP(cf,prob = prob, options = options)
print(result)
matplot(cbind(result0$x, result$x), cbind(result0$cdf, result$cdf),
        xlab = 'x', ylab = 'CDF',
        main = expression(paste('CDFs of -log(',Lambda,') under null and alternative hypothesis')))
matplot(cbind(result0$x, result$x), cbind(result0$pdf, result$pdf),
        xlab = 'x', ylab = 'CDF',
        main = expression(paste('PDFs of -log(',Lambda,') under null and alternative hypothesis')))

## EXAMPLE 4 (Compare exact vs. simulated non-central Wilk's distribution)
# p <- 10 # p - length of the sample vectors (dimensionality)
# n <- 30 # n - sample size / degrees of freedon
# q <- 5  # q - degrees of freedom due to hypothesis model
# N <- 10000  # N - number of simulation samples
# M <- replicate(q, runif(p)) # M - the (true) mean (p x q)-matrix
# delta <- eigen(M %*% t(M))$values # delta - the eigenvalues of non-centrality matrix
# L <- rep(0, N)
# for(i in 1:N) {
#         X <- replicate(n, rnorm(p))
#         E <- X %*% t(X)
#         Y <- replicate(q, rnorm(p)) + M
#         H <- Y %*% t(Y)
#         L[i] <- det((E + H) / E)
# }
# # Exact and the empirical CDF of -log(L)
# cf <- function(t) cf_LogRV_WilksLambdaNC(t, p, n, q, delta, -1)
# options <- list()
# options$xMin <- 0
# options$SixSigmaRule <- 6
# prob <- c(0.9, 0.95, 0.99)
# result <- cf2DistGP(cf, prob = prob, options = options)
# Fn <- ecdf(L)
# matplot(cbind(result$x, result$x), cbind(result$cdf, Fn(result$x)),
#         xlab = "x", ylab = "CDF / ECDF",
# title = "Exact vs. empirical CDF of the non-central distribution")

