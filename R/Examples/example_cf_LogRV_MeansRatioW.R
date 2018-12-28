## EXAMPLE 1
# CF of log MeansRatio RV with n = 5 and alpha = 7/2
n <- 5
alpha <- 7/2
t <- seq(-100, 100, length.out = 201)
plotReIm(function(t) cf_LogRV_MeansRatioW(t, n, alpha), t,
         title = 'CF of log MeansRatio RV with n = 5 and alpha = 7/2')

## EXAMPLE 2
# CF of a weighted linear combination of minus log MeansRatio RVs
rm(list=ls())
n <- c(5, 7, 10)
alpha <- list(c(7/2), c(10/2), c(3/2))
weight = list()
coef <- -1/3
t <- seq(-100, 100, length.out = 201)
plotReIm(function(t) cf_LogRV_MeansRatioW(t, n, alpha, weight, coef), t,
         title = 'CF of a weighted linear combination of minus log MeansRatio RVs')

## EXAMPLE 3
# PDF/CDF of minus log MeansRatio RV, n = 5 and alpha = 7/2, from its CF
rm(list=ls())
n <- 5
alpha <- 7/2
coef <- -1
cf <- function(t) cf_LogRV_MeansRatioW(t, n, alpha, coef = coef)
x <- seq(0, 0.6, length.out = 100)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, x, prob, options)
str(result)

## EXAMPLE 4
# PDF/CDF of minus log of the (weighted) means ratio RV, from its CF
rm(list=ls())
n <- 5
alpha <- list(c(3, 5, 7, 10, 3) / 2)
weight <- list(alpha[[1]] / sum(alpha[[1]]))
coef <- -1
cf <- function(t) cf_LogRV_MeansRatioW(t, n, alpha, weight, coef)
x <- seq(-0.15, 1.15, length.out = 200)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2^12
result <- cf2DistGP(cf, prob = prob, options = options)
str(result)

## EXAMPLE 5
# Compare the exact distribution with the Bartlett's approximation
rm(list=ls())
k <- 15 # k normal populations with unequal sample sizes
df <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3) # degrees of freedom
DF <- sum(df)
alpha <- list(df / 2)
weight <- list(alpha[[1]] / sum(alpha[[1]]))
C <- DF * log(k * prod(weight[[1]] ^ weight[[1]]))
C_B <- 1 + 1 / (3 * (k - 1)) * (sum(1 / df) - 1 / DF)
shift <- C / C_B
coef <- -DF / C_B
cf_R <- function(t) cf_LogRV_MeansRatioW(t, k, alpha, weight, coef)
cf <- function(t) exp(1i * t * shift) * cf_R(t)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
options$N <- 2^12
result <- cf2DistGP(cf, prob = prob, options = options)
str(result)
x <- result$x
matplot(cbind(x, x), cbind(result$cdf, pchisq(x, k - 1)),
        xlab = 'corrected test statistic', ylab = 'CDF',
        main = 'Exact CDF vs. the Bartlett approximation')
matplot(cbind(x, x), cbind(result$cdf, pchisq(x, k - 1)),
        xlab = 'corrected test statistic', ylab = 'CDF',
        main = 'Exact CDF vs. the Bartlett approximation',
        type = "l", lwd = 2)
print(prob)
print(result$qf)
print(qchisq(prob, k - 1))
