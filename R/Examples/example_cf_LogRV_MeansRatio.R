## EXAMPLE 1
# CF of log MeansRatio RV with n = 5 and alpha = 7/2
n <- 5
alpha <- 7/2
t <- seq(-100, 100, length.out = 201)
plotReIm(function(t) cf_LogRV_MeansRatio(t, n, alpha), t,
         title = 'CF of log MeansRatio RV with n = 5 and alpha = 7/2')

## EXAMPLE 2
# CF of a weighted linear combination of minus log MeansRatio RVs
n <- c(5, 7, 10)
alpha <- c(7, 10, 15) / 2
coef <- -c(5, 7, 10) / 22
t <- seq(-100, 100, length.out = 201)
plotReIm(function(t) cf_LogRV_MeansRatio(t, n, alpha, coef), t,
         title = 'CF of a weighted linear combination of -log MeansRatio RVs')

## EXAMPLE 3
# PDF/CDF of minus log MeansRatio RV, n = 5 and alpha = 7/2, from its CF
n <- 5
alpha <- 7/2
coef <- -1
cf <- function(t) cf_LogRV_MeansRatio(t, n, alpha, coef)
x <- seq(0, 0.6, length.out = 100)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, x, prob, options)
str(result)

## EXAMPLE 4
# Compare the exact distribution with the Bartlett's approximation
k <- 25 # number of normal populations
df <- 3  # degrees of freedom used in each of n populations
DF <- k * df
alpha <- df / 2
C_B <- (1 + 1 / (3 * (k - 1)) * (k / df - 1 / DF))
coef <- -DF / C_B
cf <- function(t) cf_LogRV_MeansRatio(t, k, alpha, coef)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, prob = prob, options = options)
str(result)
x <- result$x
matplot(cbind(x, x), cbind(result$cdf, pchisq(x, k - 1)),
        xlab = 'corrected test statistic', ylab = 'CDF',
        main = 'Exact CDF vs. the Bartlett approximation')
matplot(cbind(x, x), cbind(result$cdf, pchisq(x, k - 1)),
        xlab = 'corrected test statistic', ylab = 'CDF',
        main = 'Exact CDF vs. the Bartlett approximation', type = "l", lwd = 2)
prob
result$qf
qchisq(prob, k-1)

##  EXAMPLE 5
# Exact Critical Values for Bartlett's Test for Homogeneity of Variances
# See and compare the selected results in Glaser (1976b, Table 1)
k <- 3
df <- c(4, 5, 6, 7, 8, 9, 10, 11, 14, 19, 24, 29, 49, 99)
alpha <- df / 2
coef <- -1
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2^12
options$SixSigmaRule <- 15
options$xMin <- 0
options$isPlot <- FALSE
critW = matrix(0, length(df), length(prob))
for(i in 1:length(df)) {
        cf <- function(t) cf_LogRV_MeansRatio(t, k, alpha[i], coef)
        result <- cf2DistGP(cf, prob = prob, options = options)
        critW[i,] <- result$qf
}
critR <- exp(-critW)
cat('k = ', k, '\n')
cat('alpha = ', 1 - prob, '\n')
cat(critR)
