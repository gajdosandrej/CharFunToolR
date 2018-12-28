## EXAMPLE1 (CF of the Exponential distribution with mu = 1)
pdfFun <- function(x) exp(-x)
t <- seq(from = -20, to = 20, length.out = 2^10+1)
plotReIm(function(t) cfX_PDF(t, pdfFun), t,
         'Characteristic function of the Exponential distribution')

## EXAMPLE2 (CF of the LogNormal distribution with mu = 0, sigma = 1)
mu <- 0
sigma <- 1
pdfFun <- function(x) exp(-0.5 * ((log(x) - mu) / sigma) ^ 2) / (x * sqrt(2 * pi) * sigma)
method <- 'fit'
t <- seq(from = -20, to = 20, length.out = 2^10+1)
plotReIm(function(t) cfX_PDF(t, pdfFun, method), t,
         'Characteristic function of the LogNormal distribution')

## EXAMPLE3 (CDF/PDF of the LogNormal distribution with mu = 0, sigma = 1)
mu <- 0
sigma <- 1
pdfFun <- function(x) exp(-0.5 * ((log(x) - mu) / sigma) ^ 2) / (x * sqrt(2 * pi) * sigma)
method <- 'fit'
cf <- function(t) cfX_PDF(t, pdfFun, method)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE4 (CF of the Weibull distribution with a = 1.5, and small b<1)
a <- 1.5
b <- 0.5
pdfFun <- function(x) (x / a) ^ (b - 1) * exp(-((x / a) ^ b)) * b / a
method <- 'fit'
t <- seq(from = -20, to = 20, length.out = 2^10+1)
plotReIm(function(t) cfX_PDF(t, pdfFun, method), t,
         'Characteristic function of the Weibull distribution')

## EXAMPLE5 (CF of the Weibull distribution with a = 1.5, and large b > 1)
a <- 1.5
b <- 3.5
pdfFun <- function(x) (x / a) ^ (b - 1) * exp(-((x / a) ^ b)) * b / a
method <- 'def'
t <- seq(from = -10, to = 10, length.out = 2^10+1)
plotReIm(function(t) cfX_PDF(t, pdfFun, method), t,
         'Characteristic function of the Weibull distribution')


