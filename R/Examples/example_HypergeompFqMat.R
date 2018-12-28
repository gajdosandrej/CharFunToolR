## EXAMPLE 1
# Evaluate 2F3^9([3 4];[5 6 7];[1 2];[8,9]) & kappa<=(4,3), |kappa|<=6
a <- t(c(3, 4))
b <- t(c(5, 6, 7))
x <- t(c(1, 2))
y <- t(c(8, 9))
alpha <- 9
MAX <- 6
lam <- c(4, 3)
result <- HypergeompFqMat(a, b, x, y, alpha, MAX, lam)

## EXAMPLE 2
# CF of minus log of noncentral Wilks Lambda RV distribution
## cf_LogRV_WilksLambdaNC is using HypergeompFqMat
p <- 10
m <- 30 # elsewhere it is denoted as n (d.f. of within SS&P)
n <- 5  # elsewhere it is denoted as q (d.f. of between SS&P)
X <- c(0, 0, 0.1, 1, 20)
coef <- -1
MAX <- 25
t <- seq(-10, 10, length.out = 201)
plotReIm(function(t)
        cf_LogRV_WilksLambdaNC(t, p, m, n, X, coef, MAX),
        t,
        title = 'CF of log of noncentral Wilks Lambda RV')

## EXAMPLE 3
# PDF/CDF of minus log Wilks Lambda RV (p=10, m=30, n=5) from its CF
p <- 10
m <- 30
n <- 5
X <- c(0, 0, 0.1, 1, 20)
MAX <- 30
coef <- -1
cf0 <- function(t) cf_LogRV_WilksLambdaNC(t, p, m, n, coef = coef, MAX = MAX)
cf <- function(t) cf_LogRV_WilksLambdaNC(t, p, m, n, X, coef, MAX)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result0 <- cf2DistGP(cf0, prob = prob, options = options)
result <- cf2DistGP(cf, prob = prob, options = options)
print(result)
plot(result0$x, result0$cdf, type="l", col="orange",
     xlab = 'x', ylab = 'CDF',
     main = expression(paste('CDFs of -log(',Lambda,') under null and alternative hypothesis')))
lines(result$x, result$cdf, col="green")
plot(result0$x, result0$pdf, type="l", col="yellow",
     xlab = 'x', ylab = 'PDF',
     main = expression(paste('PDFs of -log(',Lambda,') under null and alternative hypothesis')))
lines(result$x, result$pdf, col="pink")

## EXAMPLE 4 (!!!LONG CALCULATION of Hypergeom2F1Mat!!!)
# Non-null distribution of log-transformed T statistic, W = -log(T)
# Test statistic for testing equality of 2 covariance matrices
p <- 5          # p - length of the sample vectors (dimensionality)
n1 <- 15        # n1 - sample size from the first normal poulation
n2 <- 20        # n1 - sample size from the second normal poulation
nu1 <- n1 - 1   # nu1 - degrees of freedom
nu2 <- n2 - 1   # nu1 - degrees of freedom
nu <- nu1 + nu2 # nu1 - degrees of freedom
# delta - eigenvalues of the non-centrality matrix Delta
delta <- c(1.1, 1.325, 1.55, 1.775, 2.0)
# Create the characteristic function of null distribution
c <- GammaMultiLog(nu / 2, p) - GammaMultiLog(nu1 / 2, p) - GammaMultiLog(nu2 / 2, p)
cfH0 <- function(t) {exp(c + GammaMultiLog((1 - 1i * t) * nu1 / 2, p) +
GammaMultiLog((1 - 1i * t) * nu2 / 2, p) - GammaMultiLog((1 - 1i * t) * nu / 2, p))}
# Create the characteristic function of non-null distribution
MAX <- 50
cfHA <- function(t) {cfH0(t) * prod(delta)^(-nu1 / 2) *
Hypergeom2F1Mat(nu / 2, (1 - 1i * t) * nu1 / 2, (1 - 1i * t) * nu / 2, 1 - 1 / delta, MAX)[[1]]}
# Evaluate PDF/CDF and selected quantiles of W = -log(T)
x <- seq(50, 90, length.out = 100)
prob <- c(0.9, 0.95, 0.975, 0.99, 0.999)
result <- cf2DistGP(cfHA, x, prob)
