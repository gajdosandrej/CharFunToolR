## EXAMPLE 1
# CF of log Wilks Lambda RV distribution Lambda(p,n,q)
p <- 5
n <- 10 # d.f. of within SS&P
q <- 3 # d.f. of between SS&P
t <- seq(-10, 10, length.out = 201)
plotReIm(function(t)
        cf_LogRV_WilksLambda(t, p, n, q),
        t,
        title = 'CF of log Wilks Lambda RV with Lambda(p,n,q), p=5, n=10, q=3')

## EXAMPLE 2
# CF of a weighted linear combination of minus log Wilks Lambda RVs
p <- c(5, 5, 5)
n <- c(10, 15, 20)
q <- c(3, 2, 1)
coef <- -c(10, 15, 20) / 45
t <- seq(-20, 20, length.out = 201)
plotReIm(function(t)
        cf_LogRV_WilksLambda(t, p, n, q, coef),
        t,
        title = 'CF of a weighted linear combination of -log Wilks Lambda RVs')

## EXAMPLE 3
# PDF/CDF of minus log Wilks Lambda RV (p=5, n=10, q=3) from its CF
p <- 5
n <- 10
q <- 3
coef <- -1
cf <- function(t) {cf_LogRV_WilksLambda(t, p, n, q, coef)}
x <- seq(0, 5, length.out = 100)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, x, prob, options)
str(result)

## EXAMPLE 4
# Compare the exact distribution with the Bartlett's approximation
# The Bartlett's approximation (see e.g. Wikipedia) is given by:
# ((p-q+1)/2 - n)*log(Lambda(p,n,q)) ~ chi^2_{q*p}
p <- 15
n <- 30
q <- 3
coef <- (p - q + 1) / 2 - n
cf <- function(t) {cf_LogRV_WilksLambda(t, p, n, q, coef)}
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf, prob = prob, options = options)
str(result)
x <- result$x
matplot(cbind(x, x), cbind(result$cdf, pchisq(x, q * p)),
        main = 'Exact CDF vs. the Bartlett approximation',
        xlab = expression(paste("transformed ", Lambda, "(p,n,q)")),
        ylab = 'CDF')
matplot(cbind(x, x), cbind(result$cdf, pchisq(x, q * p)),
        main = 'Exact CDF vs. the Bartlett approximation',
        xlab = expression(paste("transformed ", Lambda, "(p,n,q)")),
        ylab = 'CDF', lwd = 2, type = "l")
prob
result$qf
qchisq(prob, q * p)
