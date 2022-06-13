## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(CharFunToolR)
cf <- function(t) exp(-t ^ 2 / 2)

## ----fig.height = 4, fig.width = 6--------------------------------------------
result <- cf2DistGP(cf)
# to see the structure of the whole output you can run the following command 
# str(result)

## -----------------------------------------------------------------------------
n <- 25
p <- 0.3
t <- seq(-15, 15, length.out = 1001)

## ----fig.height = 4, fig.width = 6--------------------------------------------
plotReIm(function(t)
        cfN_Binomial(t, n, p),
        t,
        title = "CF of the Binomial distribution with n = 25, p = 0.3")

## -----------------------------------------------------------------------------
n <- 25
p <- 0.3
lambda <- 10
cfX <- function(t) cfX_Exponential(t, lambda)
t <- seq(-10, 10, length.out = 501)

## ----fig.height = 4, fig.width = 6--------------------------------------------
plotReIm(function(t)
        cfN_Binomial(t, n, p, cfX),
        t,
        title = "CF of the compound Binomial-Exponential distribution")

## -----------------------------------------------------------------------------
n <- 25
p <- 0.3
lambda <- 5
cfX <- function(t)
        cfX_Exponential(t, lambda)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)

## -----------------------------------------------------------------------------
x <- seq(0, 5, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE

## ----fig.height = 4, fig.width = 6--------------------------------------------
result <- cf2DistGP(cf, x, prob, options)
# str(result)

## ----fig.height = 4, fig.width = 6--------------------------------------------
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
plot(
        1:50,
        coef,
        xlab = "",
        ylab = "",
        type = "p",
        pch = 20,
        col = "blue",
        cex = 1,
        main = expression('Coefficients of the linear combination of GAMMA RVs')
)
lines(1:50, coef, col = "blue")

## -----------------------------------------------------------------------------
alpha <- 5 / 2
beta <- 1 / 2
t <- seq(from = -10,
         to = 10,
         length.out = 201)

## ----fig.height = 4, fig.width = 6--------------------------------------------
plotReIm(function(t)
        cf_Gamma(t, alpha, beta, coef),
        t,
        title = "Characteristic function of the linear combination of GAMMA RVs")

## -----------------------------------------------------------------------------
alpha <- 5 / 2
beta <- 1 / 2
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
cf <- function(t)
        cf_Gamma(t, alpha, beta, coef)

## -----------------------------------------------------------------------------
options <- list()
options$N <- 2 ^ 10
options$xMin <- 0

## ----fig.height = 4, fig.width = 6--------------------------------------------
result <- cf2DistGP(cf = cf, options = options)

## ----fig.height = 4, fig.width = 6--------------------------------------------
set.seed(101)
n <- 1000
data <- c(rnorm(3 * n, 5, 0.2), rt(n, 3), rchisq(n, 1))
t <- seq(-50, 50, length.out = 2 ^ 10)
plotReIm(function(t)
        cfE_Empirical(t, data),
        t,
        title = "Empirical CF - CF of the mixture of Dirac random variables")

## ----fig.height = 4, fig.width = 6--------------------------------------------
set.seed(101)
lambda <- 25
nN <- 10
Ndata <- rpois(nN, lambda)

mu <- 0.1
sigma <- 2
nX <- 1500
Xdata <- rlnorm(nX, mu, sigma)
cfX <- function(t)
        cfE_Empirical(t, Xdata)
cf  <- function(t)
        cfE_Empirical(t, Ndata, cfX)
t <- seq(-0.2, 0.2, length.out = 2 ^ 10)
plotReIm(cf, t, title = "Compound Empirical CF")

## ----fig.height = 4, fig.width = 6--------------------------------------------
x <- seq(0, 1000, length.out = 501)
prob <- c(0.9, 0.95)
options <- list()
options$N <- 2 ^ 10
options$SixSigmaRule <- 10
result <- cf2DistGP(cf, x, prob, options)

